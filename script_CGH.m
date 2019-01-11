% Script to produce a Computer Generated Hologram
% Simulation of the interference process and intensity recording on a CCD/CMOS sensor
% Reconstruction of the hologram (filter + digital refocusing)

%% Physical values (to adjust)
% spatial sampling interval on the CCD sensor, in meters
sampling = 9.3e-6;
% laser wavelength, in meters
lambda = 632.8e-9;
% distance of the object from the CCD sensor, in meters
d = 40e-2;
% reference beam tilt in x and y directions, in degrees
reference_tilt_x = 1;
reference_tilt_y = 1;
% intensity ratio between object beam and reference beam (less than or
% equal to 1 is better)
ratio = 0.95;

%% Reconstruction parameters
% reconstruction beam tilt in x and y directions, in degrees
% (should match the configuration of the reference beam)
reconstruction_tilt_x = reference_tilt_x;
reconstruction_tilt_y = reference_tilt_y;
% size of the filter used to block the zero order and twin image
filter_size_x = 200;
filter_size_y = 200;

%% Display parameters
% threshold value used for logarithmic scaling to improve contrast in the
% Fourier domain
threshold = 1e0;

%% Object slide
% replace with any greyscale picture
% transparency slide modulating the amplitude only
object_slide = double(imread('USAF.tif'));

%% Figures
% check if figures already exist from a previous run, if not create them
if isempty(findobj('Tag', 'intensity'))
    hfig_intensity = figure('Name', 'Intensity distributions', 'Tag', 'intensity');
    colormap(gray(256))
end
if isempty(findobj('Tag', 'phase'))
    hfig_phase = figure('Name', 'Phase distributions', 'Tag', 'phase');
    colormap(PhaseColormap)
end
if isempty(findobj('Tag', 'Fourier'))
    hfig_Fourier = figure('Name', 'Fourier domain', 'Tag', 'Fourier');
    colormap(gray(256))
end

%% Calculate object beam at the plane of the CCD sensor
% assume a normal plane wave passing through the transparency
% using the first Rayleigh-Sommerfeld propagation integral
object_beam = AngularSpectrumPropagation(object_slide, sampling, lambda,d);
% normalize the object beam intensity values between 0 and 1 to match the
% reference beam, then apply the intensity ratio
m1 = max(abs(object_beam(:)));
object_beam = object_beam/m1*ratio;

%% Calculate reference beam at the plane of the CCD sensor
% calculate the reference beam with same size as the object beam
size_x = size(object_slide, 2);
size_y = size(object_slide, 1);
% calculate pixel indexes, with 0 at the centre
x_vector = (1:size_x) - floor(size_x/2) - 1;
y_vector = (1:size_y) - floor(size_y/2) - 1;
% calculate spatial location for each pixel, with correct sampling interval
[x, y] = meshgrid(x_vector*sampling, y_vector*sampling);
% wave vector
k = 2*pi/lambda;
% tilt in radians
tilt_x_rad = reference_tilt_x/180*pi;
tilt_y_rad = reference_tilt_y/180*pi;
% reference beam: plane wave with constant amplitude (1)
% tilts in x and/or y lead to a phase ramp
reference_beam = exp(1i*k*(x*sin(tilt_x_rad) + y*sin(tilt_y_rad)));

%% Calculate the interference between the two beams
% intensity: square of modulus
hologram = abs(object_beam + reference_beam).^2;

%% Reconstruction beam
% tilt in radians
tilt_x_rad = reconstruction_tilt_x/180*pi;
tilt_y_rad = reconstruction_tilt_y/180*pi;
% reconstruction beam: plane wave with constant amplitude (1)
% tilts in x and/or y lead to a phase ramp
reconstruction_beam = exp(1i*k*(x*sin(tilt_x_rad) + y*sin(tilt_y_rad)));

%% Re-centering and filtering
% diffraction of the reconstruction beam on the hologram
% the effect is to re-center the virtual image in the Fourier domain, if
% the reconstruction beam matches the original reference beam
% same effect can be achieved in the Fourier domain with circshift
virtual_image = hologram.*reconstruction_beam;
virtual_image_FFT = fftshift(fft2(fftshift(virtual_image)));
% filter window
mask = zeros(size_y, size_x);
mask_vector_x = (1:filter_size_x) - floor(filter_size_x/2) - 1 + floor(size_x/2);
mask_vector_y = (1:filter_size_y) - floor(filter_size_y/2) - 1 + floor(size_y/2);
mask(mask_vector_y, mask_vector_x) = 1;
% apply the filter mask to suppress the twin image and the zero order term
virtual_image_FFT = virtual_image_FFT.*mask;
virtual_image = ifftshift(ifft2(ifftshift(virtual_image_FFT)));

%% Reconstruction of the virtual image
% using the first Rayleigh-Sommerfeld propagation integral
% minus sign: reconstruction behind the hologram for the virtual image
virtual_image_refocused = AngularSpectrumPropagation(virtual_image, sampling, lambda, -d);

%% Display the resulting intensity and phase distributions
figure(hfig_intensity)
subplot(2, 2, 1)
imagesc(abs(object_slide).^2), axis image
title({'Intensity of the object beam', 'just after the transparency slide'})
subplot(2, 2, 2)
imagesc(abs(object_beam).^2), axis image
title({'Intensity of the object beam', 'at the plane of the CCD sensor'})
subplot(2, 2, 3)
imagesc(hologram), axis image
title({'Intensity of the hologram', 'at the plane of the CCD sensor'})
subplot(2, 2, 4)
imagesc(abs(virtual_image_refocused).^2), axis image
title({'Intensity of the reconstructed image', 'at distance d'})

figure(hfig_phase)
subplot(2, 2, 1)
imagesc(angle(object_beam)), axis image
title({'Phase of the object beam alone', 'at the plane of the CCD sensor'})
subplot(2, 2, 2)
imagesc(angle(reference_beam)), axis image
title({'Phase of the reference beam alone', 'at the plane of the CCD sensor'})
subplot(2, 2, 3)
imagesc(angle(virtual_image)), axis image
title({'Phase of the filtered virtual image', 'at the plane of the CCD sensor'})
subplot(2, 2, 4)
imagesc(angle(virtual_image_refocused)), axis image
title({'Phase of the refocused virtual image', 'at distance d'})

%% Display the modulus in the Fourier domain 
object_beam_FFT = fftshift(fft2(fftshift(object_beam)));
reference_beam_FFT = fftshift(fft2(fftshift(reference_beam)));
hologram_FFT = fftshift(fft2(fftshift(hologram)));

figure(hfig_Fourier)
subplot(2, 2, 1)
imagesc(log(abs(object_beam_FFT) + threshold)), axis image
title('Object beam')
subplot(2, 2, 2)
imagesc(log(abs(reference_beam_FFT) + threshold)), axis image
title('Reference beam')
subplot(2, 2, 3)
imagesc(log(abs(hologram_FFT) + threshold)), axis image
title('Hologram')
subplot(2, 2, 4)
imagesc(log(abs(virtual_image_FFT) + threshold)), axis image
title('Virtual image filtered and re-centered')


%% this is the script to perform DRPE through simulation. 
clear,clc;close all;
%Constants
% spatial separation of points in the CCD camera,in meters
sampling = 9.3e-6;

%wavelength of light used, in metres
wavelength = 632.8e-9;

% reference tilt in degrees. in x and y
reference_tilt_x = 0;
reference_tilt_y = 0;

% intensity ratio
ratio = 0.95;
delta = pi./2;

%% reconstruction parameters
% should be the same as the reference beam
reconstructed_tilt_x = reference_tilt_x;
reconstructed_tilt_y = reference_tilt_y;
% size of the filter used to block the zero order and the virtual images
filter_size_x = 511;
filter_size_y = 511;
% size of high-pass filter used to suppress the DC term
DC_cut_size_x = 5;
DC_cut_size_y = 5;

%% display parameters 
% threshold value for logarithmic scaling to improve contrast in the
% Fourier domain
threshold = 1e0;
%% read the image into matlab
% image used here as the object changes only the amplitude
object = double(imread('USAF.tif'));
% object = rgb2gray(object);

% creating a phase-modulating object
% object_slide = object./max(object(:));
% object = exp(1i.*pi.*object_slide);

%% calculation of object plane to the plane of the 1st phase mask
% we are assuming a normal plane wave passing through the object
% distance from the object to the 1st phase mask d1 is defined, in metres
d1 = 5e-2;
%propagation to the phase mask using the Angular Spectrum method
plane_m1 = AngularSpectrumPropagation(object, sampling, wavelength, d1);

%% calculation of the effects of phase mask M1
m = size(object, 1);
n = size(object, 2);
M1 = exp(1i*pi.*(2*rand(m,n)-1));
plane_m1 = plane_m1.*M1;

%% propagation of these results to the plane of the second phase mask m2
% distance between the 1st and 2nd phase mask is d2;
d2 = 5e-2; % in metres
% propagation to the plane of the second phase mask
plane_m2 = AngularSpectrumPropagation(plane_m1, sampling, wavelength, d2);

%% calculation of the effect of the second phase mask
M2 = exp(1i.*pi.*(2*rand(m,n)-1));
plane_m2 = plane_m2.*M2;

%% propagation to the plane of the CCD camera(hologram plane)
%distance d3 is the distance between the 2nd phase mask and the CCD camera
d3 = 5e-2;
% propagation to the plane of the CCD camera
plane_holo = AngularSpectrumPropagation (plane_m2, sampling, wavelength,d3);

%% normalize the object wave 
% normalize the object beam to 0 and 1 to match the reference beam, and
% then apply the intensity ratio
max1 = max(abs(plane_holo(:)));
object_beam = plane_holo/max1*ratio;

%% calculation of the reference wave at the plane of the CCD camera
% calculate the reference beam with the same size as the object beam
size_x = size(object, 2);
size_y = size(object, 1);

% calculation of the pixel indexes, 0 at the centre
x_vector = (1:size_x) - floor(size_x/2) - 1;
y_vector = (1:size_y) - floor(size_y/2) - 1;

% calculation of the spatial location for each pixel, using the correct
% sampling interval
[x,y] = meshgrid(x_vector * sampling, y_vector* sampling);

% wave vector
k = 2*pi/ wavelength;

% tilt in radians
tilt_x_rad = reference_tilt_x *180/pi;
tilt_y_rad = reference_tilt_y *180/pi;

% the reference beam is a plane wave, with a constant amplitude of 1
% the tilt in x and y lead to a phase ramp

reference_beam = exp(1i*k*(x*sin(tilt_x_rad) + y*sin(tilt_y_rad)));

%% calculation of intensity of the 2 beams
% intensity is equal to the square of the modulus
hologram = abs(object_beam + reference_beam).^2;
hologram1 = abs(object_beam + (reference_beam.*exp(-1j.*delta))).^2;
hologram2 = abs(object_beam + (reference_beam.*exp(-2j.*delta))).^2;
hologram3 = abs(object_beam + (reference_beam.*exp(-3j.*delta))).^2;

a1 = (hologram - hologram2)./2;
a2 = (hologram1 - hologram3)./2;
hologram = a1 - 1i.*a2;
%% reconstruction of the hologram
% reconstruction beam, should be the same as the reference beam
tilt_x_rad = reconstructed_tilt_x*180/pi;
tilt_y_rad = reconstructed_tilt_y*180/pi;

% the reconstruction beam now is
reconstruction_beam = exp(1i*k*(x*sin(tilt_x_rad)+ y* sin(tilt_y_rad)));

%% Re-centering and filtering
% pass the reconstructed_beam through the hologram to diffract it. Then
% centre the virtual image in the Fourier domain
virtual_image = reconstruction_beam .* hologram;
virtual_image_FFT = fftshift(fft2(fftshift(virtual_image)));

% filter windows
mask = zeros(size_y, size_x);
mask_vector_x = (1:filter_size_x) - floor(filter_size_x/2) - 1 + floor(size_x/2);
mask_vector_y = (1:filter_size_y) - floor(filter_size_y/2) - 1 + floor(size_y/2);
mask(mask_vector_y, mask_vector_x) = 1;

mask2 = ones(size_y, size_x);
mask_vector_x = (1:DC_cut_size_x) - floor(DC_cut_size_x/2) - 1 + floor(size_x/2);
mask_vector_y = (1:DC_cut_size_y) - floor(DC_cut_size_y/2) - 1 + floor(size_y/2);
mask2(mask_vector_y, mask_vector_x) = 0;

% apply the filter to suppress the zero and real images
virtual_image_FFT = virtual_image_FFT .*mask.*mask2;
virtual_image = ifftshift(ifft2(ifftshift(virtual_image_FFT)));

%% propagate the virtual image to the object plane for reconstruction
% propagate the reconstructed beam from the hologram plane to the 2nd phase
% mask. Minus sign signifies we are reconstructing the beam behind the
% hologram
plane_holo1 = AngularSpectrumPropagation (virtual_image,sampling, wavelength, -d3);

% calculation of the effects of the second phase mask
plane_m2 = plane_holo1.*conj(M2);

% propagate to the plane of phase mask 1

plane_holo2 = AngularSpectrumPropagation(plane_m2, sampling, wavelength, -d2);

% calculation of the effects of the 1st phase mask
plane_m3 = plane_holo2.*conj(M1);
% propagation to the object plane

virtual_image_refocused = AngularSpectrumPropagation(plane_m3, sampling, wavelength, -d1);

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
%% displaying the resulting amplitude and phase distributions
figure(hfig_intensity)
subplot(2, 2, 1)
imagesc(abs(object).^2), axis image
title({'Intensity of the object beam', 'just after the transparency slide'})
subplot(2, 2, 2)
imagesc(abs(object_beam).^2), axis image
title({'Intensity of the object beam', 'at the plane of the CCD sensor'})
subplot(2, 2, 3)
imagesc(abs(hologram)), axis image
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
title({'Phase of the refocused virtual image', 'at distance d3'})

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



function output_field = RS1FFT(input_field, sampling, lambda, distance, upsample, upscale)
% First Rayleigh-Sommerfeld convolution integral
% Convolution calculated in the Fourier domain (convolution theorem)
% similar to the Kirchhoff integral, except that the source direction cosine is not included
% 
% input_field: array of (complex) optical field values, sampled at regular intervals
% sampling: sampling interval in meters (assuming same sampling in x and y directions)
% lambda: wavelength in meters
% distance: propagation distance in meters
% upsample: integer>=1, increase the number of samples in the output array by decreasing the sampling interval
% upscale: integer>=1, increase the number of samples in the output array by increasing the area

k = 2*pi/lambda;
[size_y, size_x] = size(input_field);
dmin = min([size_y, size_x]*sampling*upscale*tan(pi/2 - asin(lambda/2/sampling)));
if abs(distance)<dmin
	warning(['The minimum distance is +/-', num2str(dmin), ' m for the values provided.']);
end
output_field = zeros(size_y*upsample*upscale, size_x*upsample*upscale);

% [coord_x, coord_y] = meshgrid((0:(size_x*upscale))*sampling, (0:(size_y*upscale))*sampling);
[coord_x, coord_y] = meshgrid(((-size_x*upscale):(size_x*upscale-1))*sampling, ((-size_y*upscale):(size_y*upscale-1))*sampling);


size_conv_x = 2*size_x*upscale;
size_conv_y = 2*size_y*upscale;
input_field = fft2(input_field,size_conv_y,size_conv_x);
for m=1:upsample
    for n=1:upsample
        result = ifft2(CalculateKernel(sampling*((m-1)/upsample), sampling*((n-1))/upsample).*input_field);
        result = circshift(result,[size_y*(upscale-1)/2,size_x*(upscale-1)/2]);
        output_field(n:upsample:end,m:upsample:end) = result(1:(size_y*upscale),1:(size_x*upscale));
    end
end
output_field = -output_field*distance/2/pi*sampling^2;
% with (1i*k - 1./r) = 1i*k approximation:
% output_field = output_field*distance/1i/lambda*sampling^2;

function kernel = CalculateKernel(displacement_x, displacement_y)
	% subfunction for the calculation of the convolution kernel (impulse response)
    r = sqrt((coord_x + displacement_x).^2 + (coord_y + displacement_y).^2 + distance^2);
    kernel = (1i*k - 1./r).*exp(sign(distance).*1i.*k.*r)./r.^2;
	% with (1i*k - 1./r) = 1i*k approximation:
	% kernel = exp(sign(distance).*1i.*k.*r)./r.^2;
    kernel = fft2(fftshift(kernel));
end

end
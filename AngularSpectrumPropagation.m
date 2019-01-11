function output_field = AngularSpectrumPropagation(input_field, sampling, lambda, distance)
fmax = 1/2/sampling; %Nyquist frequency
[size_y, size_x] = size(input_field);
fx_vector = linspace(0, fmax, size_x/2+1);
fx_vector = [-fliplr(fx_vector), fx_vector(2:end-1)];
fy_vector = linspace(0, fmax, size_y/2+1);
fy_vector = [-fliplr(fy_vector), fy_vector(2:end-1)];
[fx_square, fy_square] = meshgrid(fx_vector.^2, fy_vector.^2);
phase_coefficients = exp(2i*pi*distance*sqrt(1/lambda^2 - fx_square - fy_square));
output_field = ifft2(fft2(input_field).*fftshift(phase_coefficients));
end
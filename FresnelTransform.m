function field = FresnelTransform(field, sampling, lambda, distance)
% Fresnel propagation with a single Fourier transform.
% For even lengths
% Output plane scaled in x and y by 1/(distance*lambda)

[y, x] = size(field);
[Xi, Yi] = meshgrid((-x/2:x/2-1), (-y/2:y/2-1));
Xo = Xi*lambda*distance/x/sampling; Yo = Yi*lambda*distance/y/sampling;
Xi = Xi*sampling; Yi = Yi*sampling;
field = fftshift(field.*exp(1i*pi/lambda/distance*(Xi.^2+Yi.^2)));
field = fftshift(fft2(field))*sampling^2;
field = field.*exp(1i*pi/lambda/distance*(Xo.^2+Yo.^2));
field = field*exp(2i*pi*distance/lambda)/1i/lambda/distance;
end


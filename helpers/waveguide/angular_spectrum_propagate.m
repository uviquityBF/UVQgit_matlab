function [Ez, xp, yp] = angular_spectrum_propagate(E, dx, dy, lambda_um, n_medium, dz_um, N)
% Scalar angular-spectrum propagation of a complex field by distance dz.
%
% USAGE:
%   [Ez, xp, yp] = angular_spectrum_propagate(E, dx, dy, lambda_um, n_medium, dz_um, N)
%
% E          [ny x nx] complex field on a uniform grid centred on (0,0)
% dx, dy     grid spacings (um)
% lambda_um  vacuum wavelength (um);  n_medium = refractive index of the medium
% dz_um      propagation distance (um); positive = along the propagation
%            direction (assuming the exp(i(kz - wt)) convention for E's phase)
% N          zero-padded FFT size (>= size(E)); the output lives on this N x N
%            grid (centred, same spacing) so a spreading beam is not clipped.
%
% Returns Ez on the padded grid and its axes xp, yp (um). Evanescent
% components are removed. The beam is assumed to fit inside the padded window.

    [ny, nx] = size(E);
    if N < max(nx, ny), error('pad size N=%d is smaller than the input grid (%d x %d).', N, ny, nx); end

    Ep = zeros(N, N);
    r0 = floor((N - ny) / 2);  c0 = floor((N - nx) / 2);
    Ep(r0 + (1:ny), c0 + (1:nx)) = E;

    xp = ((1:N) - (N + 1) / 2) * dx;                  % centred axes, native grid maps onto them
    yp = ((1:N) - (N + 1) / 2) * dy;
    if mod(N - nx, 2) || mod(N - ny, 2)
        warning('Odd size difference between pad and native grid: half-pixel misregistration.');
    end

    k = (0:N - 1);  k(k >= N / 2) = k(k >= N / 2) - N;
    fx = k / (N * dx);   fy = k / (N * dy);           % cycles/um
    [FX, FY] = meshgrid(fx, fy);
    arg = (n_medium / lambda_um)^2 - FX.^2 - FY.^2;
    H = exp(1i * 2 * pi * sqrt(max(arg, 0)) * dz_um) .* (arg > 0);

    Ez = ifft2(fft2(Ep) .* H);
end

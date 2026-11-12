function [a_sec, a_out] = taper_propagate(a0, C, neff_re, neff_im, L_sec, lambda_um)
% Coherent forward propagation of mode amplitudes through a staircase taper
% (eigenmode expansion, no reflections).
%
% USAGE:
%   [a_sec, a_out] = taper_propagate(a0, C, neff_re, neff_im, L_sec, lambda_um)
%
% a0        [k1 x m] input amplitudes in the modes of section 1 (m = number of
%           input cases; pass eye(k1) to get the full transmission matrix)
% C         {n-1} junction matrices; C{j} is [k(j+1) x k(j)] = <mode_new, mode_old>
% neff_re   {n} real effective indices of the kept modes of each section
% neff_im   {n} |imaginary| parts (leakage) -- amplitude decays as exp(-2*pi*Im*L/lambda)
% L_sec     [n] section lengths (um);  lambda_um vacuum wavelength (um)
%
% Returns
%   a_sec{j}  [k(j) x m] amplitudes at the END of section j (after its phase /
%             attenuation, before the junction to section j+1)
%   a_out     amplitudes after the last section (= a_sec{n})
% Power in mode k of section j is |a_sec{j}(k,:)|^2; the amplitudes are
% normalized so that sum(|a0|^2) is the power in the taper input.

    n = numel(L_sec);
    a_sec = cell(n, 1);
    a = a0;
    for j = 1:n
        beta = 2 * pi * neff_re{j}(:) / lambda_um;
        alpha = 2 * pi * abs(neff_im{j}(:)) / lambda_um;
        a = bsxfun(@times, a, exp((1i * beta - alpha) * L_sec(j)));
        a_sec{j} = a;
        if j < n, a = C{j} * a; end
    end
    a_out = a;
end

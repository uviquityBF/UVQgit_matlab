function [frac, ov, gram, Pm] = mode_overlap_fractions(x, y, E_in, E_modes, P_in_total, E_other)
% Fraction of incident power projected onto each of a set of waveguide modes.
%
% USAGE:
%   [frac, ov, gram, Pm] = mode_overlap_fractions(x, y, E_in, E_modes, P_in_total)
%
% x, y         integration-grid vectors (uniform, same length units)
% E_in         [ny x nx] complex input field on that grid
% E_modes      [ny x nx x nModes] complex mode fields on the same grid
% E_other     (optional) the OTHER transverse component of each mode, same size as
%              E_modes. When given, mode norms and the Gram matrix use the full
%              transverse vector inner product (E_modes.E_modes* + E_other.E_other*),
%              while the overlap with the (single-polarization) input still uses
%              only E_modes. This avoids counting nearly-degenerate hybrid modes
%              twice and keeps modes of the other polarization from being
%              normalized up from numerical noise.
% P_in_total   total incident power (integral of |E_in|^2 over the FULL input,
%              not just the integration window). Power outside the window is
%              therefore counted as unaccounted, not silently discarded.
%
% Returns
%   frac   [nModes x 1]  |<in,m>|^2 / (P_in_total * <m,m>)
%   ov     [nModes x 1]  complex overlap <in,m>
%   gram   [nModes x nModes] |<mi,mj>| / sqrt(<mi,mi><mj,mj>)  (identity if
%          the set is orthogonal -- use as a sanity check)
%   Pm     [nModes x 1]  mode norms <m,m> (vector norm if E_other given)
%
% For an orthogonal mode set, sum(frac) <= 1 (Bessel); the deficit is power in
% modes not supplied (radiation, cladding, other polarization, window clip).

    integ = @(F) trapz(y, trapz(x, F, 2));
    n = size(E_modes, 3);
    if nargin < 6 || isempty(E_other), E_other = zeros(size(E_modes)); end

    ov = zeros(n, 1);  Pm = zeros(n, 1);
    for k = 1:n
        Ek = E_modes(:, :, k);
        ov(k) = integ(E_in .* conj(Ek));
        Pm(k) = real(integ(abs(Ek).^2 + abs(E_other(:, :, k)).^2));
    end
    frac = abs(ov).^2 ./ (P_in_total * Pm);

    gram = zeros(n, n);
    for i = 1:n
        for j = i:n
            g = abs(integ(E_modes(:, :, i) .* conj(E_modes(:, :, j)) + ...
                          E_other(:, :, i) .* conj(E_other(:, :, j)))) / sqrt(Pm(i) * Pm(j));
            gram(i, j) = g;  gram(j, i) = g;
        end
    end
end

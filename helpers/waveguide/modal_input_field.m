function E_in = modal_input_field(P, shift_um, dz_um, extra_phase_fn)
% Complex input field sampled on the integration grid of a modal problem P
% (from load_modal_problem), after an optional defocus and lateral shift.
%
% USAGE:
%   E_in = modal_input_field(P, [dx dy], dz_um, extra_phase_fn)
%
% shift_um        [dx dy] total shift of the beam on the chip grid (um)
% dz_um           propagate the POP field by this distance first (0 = as exported)
% extra_phase_fn  optional @(X,Y) radians, X,Y relative to the shifted beam centre

    if nargin < 3 || isempty(dz_um), dz_um = 0; end
    if nargin < 4, extra_phase_fn = []; end

    if dz_um == 0
        Es = P.E_native;  xs = P.zX;  ys = P.zY;
    else
        dxz = abs(P.zX(2) - P.zX(1));  dyz = abs(P.zY(2) - P.zY(1));
        [Es, xs, ys] = angular_spectrum_propagate(P.E_native, dxz, dyz, ...
                          P.lambda_um, P.n_medium, dz_um, P.pad_n);
    end
    Er = interp2(xs, ys, real(Es), P.IX - shift_um(1), P.IY - shift_um(2), 'linear', 0);
    Ei = interp2(xs, ys, imag(Es), P.IX - shift_um(1), P.IY - shift_um(2), 'linear', 0);
    E_in = Er + 1i * Ei;

    if ~isempty(extra_phase_fn)
        E_in = E_in .* exp(1i * extra_phase_fn(P.IX - (P.cx_in + shift_um(1)), ...
                                               P.IY - (P.cy_in + shift_um(2))));
    end
end

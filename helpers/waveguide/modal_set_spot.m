function P = modal_set_spot(P, spot)
% Replace the input spot inside an already-loaded modal problem (the COMSOL modes,
% which are expensive to load, are kept). Use for sweeps over spot parameters.
%
% USAGE:
%   P2 = modal_set_spot(P, make_gaussian_spot(cfg));

    P.zX = spot.zX;  P.zY = spot.zY;  P.E_native = spot.E_native;
    P.P_in_total = spot.P_in_total;  P.cx_in = spot.cx_in;  P.cy_in = spot.cy_in;
    P.lambda_um = spot.lambda_um;  P.n_medium = spot.n_medium;  P.fileI = spot.name;
    P.align_shift = P.mode_centroid - [spot.cx_in, spot.cy_in];
end

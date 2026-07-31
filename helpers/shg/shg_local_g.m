function g = shg_local_g(w_nm, h_nm, geom)
% Nonlinear coupling coefficient g at a specific local waveguide geometry.
%
% w_nm, h_nm : local width and height [nm]
% geom       : struct with fields lam_nm, d33_pmV, n_pump, n_shg, overlap_eta
%              (When geom.overlap_table is present, overlap_eta is interpolated
%               from it via shg_interp_overlap(w_nm, h_nm, geom.overlap_table).)
%
% This is the same physics as shg_g_coefficient, but accepts explicit geometry
% at a point so it can be called per-step for z-dependent g(z) simulation.
    if isfield(geom, 'g_override') && ~isempty(geom.g_override)
        g = geom.g_override;
        return;
    end

    eps0   = 8.854e-12;
    c      = 2.998e8;
    lam1   = geom.lam_nm * 1e-9;
    omega1 = 2 * pi * c / lam1;
    deff   = geom.d33_pmV * 1e-12;
    n1     = geom.n_pump;
    n2     = geom.n_shg;
    A_eff  = (w_nm * 1e-9) * (h_nm * 1e-9);

    if isfield(geom, 'overlap_table') && ~isempty(geom.overlap_table)
        ol_eta = shg_interp_overlap(w_nm, h_nm, geom.overlap_table);
    else
        ol_eta = geom.overlap_eta;
    end

    ol_factor = ol_eta / sqrt(A_eff);
    yeta      = 2 * omega1^2 / (eps0 * c^3 * n1^2 * n2) * deff^2 * ol_factor^2;
    g         = sqrt(yeta) * 1e-6;
end

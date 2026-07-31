function g = shg_g_coefficient(cfg)
% Nonlinear coupling coefficient g (Wei's formula) from a case config struct.
% cfg.phys must have: lam_nm, d33_pmV, n_pump, n_shg, width_um, height_um, overlap_eta.
    eps0   = 8.854e-12;
    c      = 2.998e8;
    lam1   = cfg.phys.lam_nm * 1e-9;
    omega1 = 2 * pi * c / lam1;
    deff   = cfg.phys.d33_pmV * 1e-12;
    n1     = cfg.phys.n_pump;
    n2     = cfg.phys.n_shg;
    A_eff  = (cfg.phys.width_um * cfg.phys.height_um) * 1e-12;

    ol_factor = cfg.phys.overlap_eta / sqrt(A_eff);
    yeta      = 2 * omega1^2 / (eps0 * c^3 * n1^2 * n2) * deff^2 * ol_factor^2;
    g         = sqrt(yeta) * 1e-6;
end

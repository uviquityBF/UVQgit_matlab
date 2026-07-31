function [Pg, Ps] = shg_run_core_simulation(cfg)
% Dispatches to the single-segment or two-segment RK4 engine based on whether
% cfg.phys.L1_um is present and non-empty (set by shg_resolve_taper_params for
% cases that have a two-segment taper struct).
%
% g is now computed as a function of z via shg_local_g, using the local waveguide
% width and height at each RK4 step. For cases without a taper, width and height
% are constant (g is effectively uniform). overlap_eta is still a constant scalar
% until cfg.phys.overlap_table is provided (see shg_local_g).
    Pp  = cfg.phys.Pp_peak_mW * 1e-3;
    a0  = db_to_np_per_um(cfg.phys.a0_dBcm);
    as  = db_to_np_per_um(cfg.phys.a3_scat_dBcm);
    aa  = db_to_np_per_um(cfg.phys.a3_abs_dBcm);
    a3t = as + aa;

    a0_grad  = db_to_np_per_um(cfg.phys.a0_grad_dBcm_mm) * 1e-3;
    a3t_grad = db_to_np_per_um(cfg.phys.a3_grad_dBcm_mm) * 1e-3;
    if a3t > 0
        a3s_grad = a3t_grad * (as / a3t);
    else
        a3s_grad = 0;
    end

    % Default geometry for g(z): uniform waveguide at the case-level dimensions.
    w_nom_nm  = cfg.phys.width_um  * 1000;
    h0_nm     = cfg.phys.height_um * 1000;
    dh_per_um = 0;

    if isfield(cfg.phys, 'L1_um') && ~isempty(cfg.phys.L1_um)
        L1 = cfg.phys.L1_um;

        % Width geometry: use promoted fields from shg_resolve_taper_params
        % (w_start_nm, w_L1_nm, w_end_nm) when available; fall back to constant
        % width for legacy rate-spec cases that don't set w_start_nm.
        if isfield(cfg.phys, 'w_start_nm') && ~isempty(cfg.phys.w_start_nm)
            w_seg1_start = cfg.phys.w_start_nm;
        else
            w_seg1_start = w_nom_nm;   % legacy: constant width throughout
        end
        w_seg1_end   = cfg.phys.w_L1_nm;
        w_seg2_start = cfg.phys.w_L1_nm;

        if isfield(cfg.phys, 'dh_per_um') && ~isempty(cfg.phys.dh_per_um)
            dh_per_um = cfg.phys.dh_per_um;
        end

        % Seg2 dk_grad and exit width: width-endpoint spec computes per current L;
        % legacy rate-spec has dk_grad_2 pre-stored and uses constant width.
        if isfield(cfg.phys, 'w_end_nm') && ~isempty(cfg.phys.w_end_nm)
            L2 = cfg.phys.L_um - L1;
            if L2 > 0
                dw_mm_2 = (cfg.phys.w_end_nm - cfg.phys.w_L1_nm) / (L2 * 1e-3);
            else
                dw_mm_2 = 0;
            end
            dk_grad_2  = cfg.phys.dk_grad_h + cfg.phys.C_scale_dw * dw_mm_2;
            w_seg2_end = cfg.phys.w_end_nm;
        else
            dk_grad_2  = cfg.phys.dk_grad_2;
            w_seg2_end = w_nom_nm;   % legacy: unknown endpoint, hold constant
        end

        geom1 = make_geom(w_seg1_start, w_seg1_end,   h0_nm, dh_per_um, cfg.phys);
        geom2 = make_geom(w_seg2_start, w_seg2_end,   h0_nm, dh_per_um, cfg.phys);

        [Pg, Ps] = shg_rk4_engine_segmented( ...
            cfg.phys.L_um, L1, cfg.sim.dz, geom1, geom2, ...
            a0, a0_grad, a3t, a3t_grad, as, a3s_grad, ...
            cfg.phys.dk_center, cfg.phys.dk_grad, dk_grad_2, Pp);
    else
        geom = make_geom(w_nom_nm, w_nom_nm, h0_nm, dh_per_um, cfg.phys);

        [Pg, Ps] = shg_rk4_engine( ...
            cfg.phys.L_um, cfg.sim.dz, geom, ...
            a0, a0_grad, a3t, a3t_grad, as, a3s_grad, ...
            cfg.phys.dk_center, cfg.phys.dk_grad, Pp);
    end
end


function geom = make_geom(w_start_nm, w_end_nm, h0_nm, dh_per_um, phys)
% Build a geometry profile struct for one waveguide segment.
    geom.w_start_nm  = w_start_nm;
    geom.w_end_nm    = w_end_nm;
    geom.h0_nm       = h0_nm;
    geom.dh_per_um   = dh_per_um;
    geom.lam_nm      = phys.lam_nm;
    geom.d33_pmV     = phys.d33_pmV;
    geom.n_pump      = phys.n_pump;
    geom.n_shg       = phys.n_shg;
    geom.overlap_eta = phys.overlap_eta;
    % Future: if isfield(phys, 'overlap_table'), geom.overlap_table = phys.overlap_table; end
end

function [Pg, Ps] = shg_run_core_simulation(cfg)
% Dispatches to the single-segment or two-segment RK4 engine based on whether
% cfg.phys.L1_um is present and non-empty (set by shg_resolve_taper_params for
% cases that have a two-segment taper struct).
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

    if isfield(cfg.phys, 'L1_um') && ~isempty(cfg.phys.L1_um)
        L1 = cfg.phys.L1_um;

        % Seg2 dk_grad: width-endpoint spec computes per current L; old dw_per_mm spec is pre-stored.
        if isfield(cfg.phys, 'w_end_nm') && ~isempty(cfg.phys.w_end_nm)
            L2 = cfg.phys.L_um - L1;
            if L2 > 0
                dw_mm_2 = (cfg.phys.w_end_nm - cfg.phys.w_L1_nm) / (L2 * 1e-3);
            else
                dw_mm_2 = 0;
            end
            dk_grad_2 = cfg.phys.dk_grad_h + cfg.phys.C_scale_dw * dw_mm_2;
        else
            dk_grad_2 = cfg.phys.dk_grad_2;
        end

        [Pg, Ps] = shg_rk4_engine_segmented( ...
            cfg.phys.L_um, L1, cfg.sim.dz, cfg.phys.g, ...
            a0, a0_grad, a3t, a3t_grad, as, a3s_grad, ...
            cfg.phys.dk_center, cfg.phys.dk_grad, dk_grad_2, Pp);
    else
        [Pg, Ps] = shg_rk4_engine( ...
            cfg.phys.L_um, cfg.sim.dz, cfg.phys.g, ...
            a0, a0_grad, a3t, a3t_grad, as, a3s_grad, ...
            cfg.phys.dk_center, cfg.phys.dk_grad, Pp);
    end
end

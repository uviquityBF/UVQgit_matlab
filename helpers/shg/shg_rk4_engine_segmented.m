function [Pg_end, Ps_end] = shg_rk4_engine_segmented(L, L1, dz, geom1, geom2, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, dk_center_1, dk_grad_1, dk_grad_2, Pp)
% Two-segment RK4 engine.
%
% geom1 : geometry profile for segment 1 (z_abs 0 -> L1).
%           geom1.w_start_nm = width at waveguide input.
%           geom1.w_end_nm   = width at segment boundary (L1).
% geom2 : geometry profile for segment 2 (z_abs L1 -> L).
%           geom2.w_start_nm = width at segment boundary (= geom1.w_end_nm).
%           geom2.w_end_nm   = width at waveguide output.
% Both geom structs share the same h0_nm and dh_per_um (height taper is global).
%
% If L <= L1, only segment 1 runs; its exit width is interpolated to L/L1
% of the full seg1 taper so g(z) remains physically correct.
%
% dk continuity: dk_center_2 = dk_center_1 + dk_grad_1 * L1
    Y0 = [sqrt(Pp); 0; 0; 0; 0];

    if L <= L1
        % Waveguide shorter than seg1 — truncate taper to actual exit width.
        geom1_run = geom1;
        if L1 > 0
            geom1_run.w_end_nm = geom1.w_start_nm + ...
                (geom1.w_end_nm - geom1.w_start_nm) * (L / L1);
        end
        [Y_end, ~] = shg_rk4_segment(L, dz, geom1_run, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, ...
                                      dk_center_1, dk_grad_1, Y0, 0.0, 0.0);
    else
        [Y1, phase1] = shg_rk4_segment(L1, dz, geom1, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, ...
                                        dk_center_1, dk_grad_1, Y0, 0.0, 0.0);
        dk_center_2 = dk_center_1 + dk_grad_1 * L1;
        L2 = L - L1;
        [Y_end, ~] = shg_rk4_segment(L2, dz, geom2, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, ...
                                      dk_center_2, dk_grad_2, Y1, L1, phase1);
    end

    Pg_end = Y_end(3)^2 + Y_end(4)^2;
    Ps_end = Y_end(5);
end

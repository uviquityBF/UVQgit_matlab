function [Pg_end, Ps_end] = shg_rk4_engine_segmented(L, L1, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, dk_center_1, dk_grad_1, dk_grad_2, Pp)
% Two-segment RK4 engine.
%
% seg1: z_local 0 -> L1,   z_abs 0 -> L1,   phase starts at 0
% seg2: z_local 0 -> L-L1, z_abs L1 -> L,   phase starts at phase_end_seg1
%
% dk_center_2 = dk_center_1 + dk_grad_1 * L1  (dk is continuous at boundary)
%
% If L <= L1, only seg1 runs (waveguide shorter than the first segment).
    Y0 = [sqrt(Pp); 0; 0; 0; 0];

    if L <= L1
        [Y_end, ~] = shg_rk4_segment(L, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, ...
                                      dk_center_1, dk_grad_1, Y0, 0.0, 0.0);
    else
        [Y1, phase1] = shg_rk4_segment(L1, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, ...
                                        dk_center_1, dk_grad_1, Y0, 0.0, 0.0);
        dk_center_2 = dk_center_1 + dk_grad_1 * L1;
        L2 = L - L1;
        [Y_end, ~] = shg_rk4_segment(L2, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, ...
                                      dk_center_2, dk_grad_2, Y1, L1, phase1);
    end

    Pg_end = Y_end(3)^2 + Y_end(4)^2;
    Ps_end = Y_end(5);
end

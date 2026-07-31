function [Pg_end, Ps_end] = shg_rk4_engine(L, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, dk_center, dk_grad, Pp)
% Single-segment RK4 engine (backward-compatible wrapper around shg_rk4_segment).
    Y0 = [sqrt(Pp); 0; 0; 0; 0];
    [Y_end, ~] = shg_rk4_segment(L, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, ...
                                  dk_center, dk_grad, Y0, 0.0, 0.0);
    Pg_end = Y_end(3)^2 + Y_end(4)^2;
    Ps_end = Y_end(5);
end

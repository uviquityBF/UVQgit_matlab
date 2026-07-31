function [Y_end, phase_end] = shg_rk4_segment(L_seg, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, dk_center, dk_grad, Y0, z_offset, phase_offset)
% RK4 integration over one piecewise section of the waveguide.
%
% Returns Y_end (state at end of segment) and phase_end (total accumulated
% phase at end of segment, to be passed as phase_offset to the next segment).
    z_vec = 0:dz:L_seg;
    Y     = Y0;

    for i = 1:(length(z_vec) - 1)
        z_loc = z_vec(i);
        z_abs = z_loc + z_offset;
        k1 = shg_cme_derivatives(z_loc,       z_abs,       Y,           a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset);
        k2 = shg_cme_derivatives(z_loc+dz/2, z_abs+dz/2, Y+k1*dz/2,   a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset);
        k3 = shg_cme_derivatives(z_loc+dz/2, z_abs+dz/2, Y+k2*dz/2,   a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset);
        k4 = shg_cme_derivatives(z_loc+dz,   z_abs+dz,   Y+k3*dz,     a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset);
        Y  = Y + (dz/6) * (k1 + 2*k2 + 2*k3 + k4);
    end

    Y_end     = Y;
    phase_end = phase_offset + dk_center * L_seg + 0.5 * dk_grad * L_seg^2;
end

function [Y_end, phase_end] = shg_rk4_segment(L_seg, dz, geom, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, dk_center, dk_grad, Y0, z_offset, phase_offset)
% RK4 integration over one piecewise section of the waveguide.
%
% L_seg        : length of this segment [um]
% dz           : step size [um]
% geom         : geometry profile for this segment (see shg_local_g).
%                  geom.w_start_nm  -- width at z_loc = 0      [nm]
%                  geom.w_end_nm    -- width at z_loc = L_seg  [nm]
%                  geom.h0_nm       -- height at z_abs = 0     [nm]
%                  geom.dh_per_um   -- height taper rate        [nm/um]
%                  + lam_nm, d33_pmV, n_pump, n_shg, overlap_eta
%                Width is interpolated linearly across the segment.
%                Height uses absolute z (continuous across segment boundaries).
% z_offset     : z_abs at the start of this segment [um]
% phase_offset : phase accumulated in all prior segments [rad]
%
% g is computed once per step at z_loc = z_vec(i) — a step-constant
% approximation that is accurate when dz is small relative to the taper length.
    z_vec = 0:dz:L_seg;
    Y     = Y0;

    for i = 1:(length(z_vec) - 1)
        z_loc = z_vec(i);
        z_abs = z_loc + z_offset;

        % Local geometry: linear width taper within segment; global height taper.
        if L_seg > 0
            w_nm = geom.w_start_nm + (geom.w_end_nm - geom.w_start_nm) * (z_loc / L_seg);
        else
            w_nm = geom.w_start_nm;
        end
        h_nm = geom.h0_nm + geom.dh_per_um * z_abs;
        g    = shg_local_g(w_nm, h_nm, geom);

        k1 = shg_cme_derivatives(z_loc,       z_abs,       Y,           a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset);
        k2 = shg_cme_derivatives(z_loc+dz/2, z_abs+dz/2, Y+k1*dz/2,   a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset);
        k3 = shg_cme_derivatives(z_loc+dz/2, z_abs+dz/2, Y+k2*dz/2,   a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset);
        k4 = shg_cme_derivatives(z_loc+dz,   z_abs+dz,   Y+k3*dz,     a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset);
        Y  = Y + (dz/6) * (k1 + 2*k2 + 2*k3 + k4);
    end

    Y_end     = Y;
    phase_end = phase_offset + dk_center * L_seg + 0.5 * dk_grad * L_seg^2;
end

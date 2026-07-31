function dY = shg_cme_derivatives(z_loc, z_abs, Y, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset)
% Coupled-mode-equation derivatives for one RK4 sub-step, within one piecewise segment.
%
% z_loc        : position local to this segment [um]  -- drives phase
% z_abs        : absolute position from waveguide input [um]  -- drives loss gradient
% phase_offset : phase accumulated in all prior segments [rad]
%
% Phase: Phi = phase_offset + dk_center*z_loc + 0.5*dk_grad*z_loc^2
% Loss:  a(z) = a_center + a_grad * z_abs   (continuous across segments)
    A1 = Y(1) + 1i*Y(2);
    A3 = Y(3) + 1i*Y(4);

    phase_accum = phase_offset + dk_center * z_loc + 0.5 * dk_grad * z_loc^2;

    a0_z  = a0  + a0_grad  * z_abs;
    a3t_z = a3t + a3t_grad * z_abs;
    a3s_z = a3s + a3s_grad * z_abs;

    dA1 = -0.5*a0_z*A1;
    dA3 = -0.5*a3t_z*A3 + 1i*g*(A1^2)*exp(-1i*phase_accum);
    dPs = a3s_z * abs(A3)^2;

    dY = [real(dA1); imag(dA1); real(dA3); imag(dA3); dPs];
end

function SHG_Efficiency_w_Scatter()
    % SHG_Master_Consolidated: RK4 Solver for TM00-TM04 SHG in AlN
    % Handles: Scattering vs. Guided power, spectral broadening, and absorption.
    
    close all; clear;

    %% 1. CONFIGURATION / INPUTS
    % Setup speed toggle
    cfg.sim.use_linewidth           = true;     % [Logic] set to false for instant runs
    cfg.sim.dz                      = 0.5;      % [um] RK4 step size
    cfg.sim.N_overlap_steps         = 51;       % Number of wavelengths for spectral integral

    % Pump Parameters
    cfg.pump.lam_nm                 = 450;
    cfg.pump.P_mW                   = 100;
    cfg.pump.linewidth_nm           = 0.15;     % FWHM of pump
    cfg.pump.loss_dBcm              = 5;

    % SHG Parameters (TM04 Mode)
    cfg.shg.loss_scat_dBcm          = 200;      % Your "Signal" (Extraction)
    cfg.shg.loss_abs_dBcm           = 20;       % Parasitic absorption
    cfg.shg.d_eff                   = 4.7;      % [pm/V]
    cfg.shg.overlap_frac            = 0.03;     % Modal overlap factor
    
    % Waveguide Geometry
    cfg.geo.L_um                    = 2000;
    cfg.geo.width_um                = 0.3365;
    cfg.geo.height_um               = 0.381;    
    cfg.geo.dispersion_slope        = 0.0089;   % [1/nm] difference in n_eff slopes

    % Sweep Parameters (for Bandwidth Study)
    detuning_list = linspace(-1.0, 1.0, 61);    % [nm] Detuning from phase match

    %% 2. PRE-CALC & CONVERSIONS
    Pp = cfg.pump.P_mW * 1e-3;
    a0 = (cfg.pump.loss_dBcm / 4.3429) * 1e-4;
    a3_scat = (cfg.shg.loss_scat_dBcm / 4.3429) * 1e-4;
    a3_abs = (cfg.shg.loss_abs_dBcm / 4.3429) * 1e-4;
    a3_total = a3_scat + a3_abs;
    
    % Calculate g (Simplified nonlinear coupling)
    % Consistent with your Version 3 logic
    w2 = 2*pi * 3e8 / ( (cfg.pump.lam_nm/2) * 1e-9);
    prefactor = 2 * w2^2 * (cfg.shg.d_eff * 1e-12)^2 / (8.85e-12 * 2.15^2 * 2.6 * 3e8^3);
    A_eff = cfg.geo.width_um * cfg.geo.height_um;
    g_coeff = sqrt( cfg.shg.overlap_frac^2 / A_eff * prefactor );

    %% 3. BANDWIDTH STUDY (Loop over Detuning)
    P_guided_detune = zeros(size(detuning_list));
    P_scat_detune = zeros(size(detuning_list));

    fprintf('Starting sweep over %d detuning points...\n', length(detuning_list));
    
    for k = 1:length(detuning_list)
        dk_center = (4*pi / (cfg.pump.lam_nm * 1e-7)) * (detuning_list(k) * cfg.geo.dispersion_slope);
        dk_um = dk_center * 1e-4;

        if cfg.sim.use_linewidth
            % --- Case A: Spectral Integration (Linewidth) ---
            HWHM = cfg.pump.linewidth_nm / 2;
            % Shifted spectrum to follow the detuning
            wl_spectrum = linspace(-3*HWHM, 3*HWHM, cfg.sim.N_overlap_steps);
            weights = exp(-0.5 * (wl_spectrum / HWHM).^2);
            weights = weights / sum(weights);
            
            p_g_temp = 0; p_s_temp = 0;
            for w_idx = 1:length(wl_spectrum)
                % Add the offset from the pump's internal spectrum
                dk_offset = (4*pi / (cfg.pump.lam_nm * 1e-7)) * (wl_spectrum(w_idx) * cfg.geo.dispersion_slope) * 1e-4;
                [pg, ps] = run_rk4(cfg, g_coeff, a0, a3_total, a3_scat, dk_um + dk_offset, Pp);
                p_g_temp = p_g_temp + pg * weights(w_idx);
                p_s_temp = p_s_temp + ps * weights(w_idx);
            end
            P_guided_detune(k) = p_g_temp;
            P_scat_detune(k) = p_s_temp;
        else
            % --- Case B: Monochromatic (Instant) ---
            [P_guided_detune(k), P_scat_detune(k)] = run_rk4(cfg, g_coeff, a0, a3_total, a3_scat, dk_um, Pp);
        end
    end

    %% 4. ANALYZE & PLOT
    % Extraction Ratio FOM
    Ratio = P_scat_detune ./ (P_guided_detune + 1e-15);
    
    % Find FWHM Broadening
    fwhm_g = calculate_fwhm(detuning_list, P_guided_detune);
    fwhm_s = calculate_fwhm(detuning_list, P_scat_detune);

    figure('Color','w','Name','SHG Power Studies');
    subplot(2,1,1);
    plot(detuning_list, P_guided_detune*1e6, 'b', 'LineWidth', 2); hold on;
    plot(detuning_list, P_scat_detune*1e6, 'r--', 'LineWidth', 2);
    grid on; xlabel('Pump Detuning [nm]'); ylabel('Power [\muW]');
    legend('Guided SHG', 'Scattered SHG');
    title(sprintf('SHG Power vs Detuning (FWHM_g: %.3f nm, FWHM_s: %.3f nm)', fwhm_g, fwhm_s));

    subplot(2,1,2);
    plot(detuning_list, Ratio, 'k', 'LineWidth', 2);
    grid on; xlabel('Pump Detuning [nm]'); ylabel('Ratio (Scat/Guided)');
    title('Extraction Ratio vs Detuning');
    
    fprintf('FWHM Guided: %.3f nm\n', fwhm_g);
    fprintf('FWHM Scattered: %.3f nm\n', fwhm_s);
    fprintf('Bandwidth Broadening Factor: %.2f\n', fwhm_s / fwhm_g);
end

%% HELPER: RK4 Inner Engine
function [Pg_end, Ps_end] = run_rk4(cfg, g, a0, a3_total, a3_scat, dk, Pp)
    z_vec = 0:cfg.sim.dz:cfg.geo.L_um;
    % Y = [A1_r, A1_i, A3_r, A3_i, P_scat, Phase]
    Y = [sqrt(Pp); 0; 0; 0; 0; 0];
    dz = cfg.sim.dz;
    
    for i = 1:(length(z_vec)-1)
        k1 = derivatives(Y, a0, a3_total, a3_scat, g, dk);
        k2 = derivatives(Y + k1*dz/2, a0, a3_total, a3_scat, g, dk);
        k3 = derivatives(Y + k2*dz/2, a0, a3_total, a3_scat, g, dk);
        k4 = derivatives(Y + k3*dz, a0, a3_total, a3_scat, g, dk);
        Y = Y + (dz/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
    Pg_end = Y(3)^2 + Y(4)^2;
    Ps_end = Y(5);
end

function dY = derivatives(Y, a0, a3_total, a3_scat, g, dk)
    A1 = Y(1) + 1i*Y(2);
    A3 = Y(3) + 1i*Y(4);
    phase = Y(6);
    dA1 = -0.5 * a0 * A1;
    dA3 = -0.5 * a3_total * A3 + 1i * g * (A1^2) * exp(-1i * phase);
    dPscat = a3_scat * (abs(A3)^2);
    dPhase = dk;
    dY = [real(dA1); imag(dA1); real(dA3); imag(dA3); dPscat; dPhase];
end

%% HELPER: FWHM Calculation
function lw = calculate_fwhm(x, y)
    y = y - min(y); % baseline
    half_max = max(y) / 2;
    idx = find(y >= half_max);
    if length(idx) < 2
        lw = 0;
    else
        lw = x(idx(end)) - x(idx(1));
    end
end
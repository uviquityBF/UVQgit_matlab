function SHG_Efficiency_w_Scatter()
    % SHG_Master_Consolidated: RK4 Solver for TM00-TM04 SHG in AlN
    % Handles: Scattering vs. Guided power, spectral broadening, and absorption.
    %
    % NOT a length/loss sweep like SHG_Design_Suite.m -- this sweeps pump
    % wavelength DETUNING at one fixed geometry/length, to compare the
    % phase-matching bandwidth (FWHM) of guided vs. scattered SHG. Reports
    % the FWHM broadening factor as the headline result.
    %
    % Shares SHG_Design_Suite.m's RK4 engine (helpers/shg/shg_rk4_engine.m),
    % but its own g-coefficient formula uses less-precise physical constants
    % (c=3e8, eps0=8.85e-12 vs. 2.998e8/8.854e-12) -- left as-is during
    % cleanup rather than silently changed, so its numbers won't exactly
    % match SHG_Design_Suite.m's at identical nominal inputs.

    close all; clear;
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\shg');

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
    a0 = db_to_np_per_um(cfg.pump.loss_dBcm);
    a3_scat = db_to_np_per_um(cfg.shg.loss_scat_dBcm);
    a3_abs = db_to_np_per_um(cfg.shg.loss_abs_dBcm);
    a3_total = a3_scat + a3_abs;
    
    % Calculate g (Simplified nonlinear coupling)
    % Consistent with your Version 3 logic
    w2 = 2*pi * 3e8 / ( (cfg.pump.lam_nm/2) * 1e-9);
    prefactor = 2 * w2^2 * (cfg.shg.d_eff * 1e-12)^2 / (8.85e-12 * 2.15^2 * 2.6 * 3e8^3);
    A_eff = cfg.geo.width_um * cfg.geo.height_um;
    g_coeff = sqrt( cfg.shg.overlap_frac^2 / A_eff * prefactor );

    % Build geom struct for shg_rk4_engine (uniform waveguide, g pre-computed above).
    % g_override bypasses shg_local_g so this script's g formula is preserved exactly.
    geom = struct('w_start_nm', cfg.geo.width_um*1000, 'w_end_nm', cfg.geo.width_um*1000, ...
                  'h0_nm', cfg.geo.height_um*1000, 'dh_per_um', 0, 'g_override', g_coeff);

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
                [pg, ps] = shg_rk4_engine(cfg.geo.L_um, cfg.sim.dz, geom, ...
                    a0, 0, a3_total, 0, a3_scat, 0, dk_um + dk_offset, 0, Pp);
                p_g_temp = p_g_temp + pg * weights(w_idx);
                p_s_temp = p_s_temp + ps * weights(w_idx);
            end
            P_guided_detune(k) = p_g_temp;
            P_scat_detune(k) = p_s_temp;
        else
            % --- Case B: Monochromatic (Instant) ---
            [P_guided_detune(k), P_scat_detune(k)] = shg_rk4_engine(cfg.geo.L_um, cfg.sim.dz, geom, ...
                a0, 0, a3_total, 0, a3_scat, 0, dk_um, 0, Pp);
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
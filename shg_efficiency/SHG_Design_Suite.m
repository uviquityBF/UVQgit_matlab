function SHG_Design_Suite()
% SHG Design Suite v3.3
%
% Changes from v3.2:
%   - Two-segment piecewise width taper.
%     Add L1_um to the 'taper' sub-struct to activate; segment 2 gets the
%     remainder of the full waveguide length.  Each segment has an independent
%     width taper rate; the height taper (dh_per_mm) applies globally.
%
%     Taper struct fields for two-segment mode:
%       cases(k).taper.L1_um        = 800;    % segment 1 fixed length [um]
%       cases(k).taper.w_start_nm   = 300;    % width at waveguide input [nm]
%       cases(k).taper.w_L1_nm      = 310;    % width at seg1/seg2 boundary [nm]
%       cases(k).taper.w_end_nm     = 310;    % width at waveguide output [nm]
%       cases(k).taper.dlam_PM_dw   = -0.012; % shared width PM sensitivity [nm/nm]
%       cases(k).taper.dh_per_mm    = 3.0;    % height taper (global) [nm/mm]
%       cases(k).taper.dlam_PM_dh   = 0.22;   % height PM sensitivity [nm/nm]
%       cases(k).taper.lam_PM_0_nm  = 450.0;  % PM wavelength at z=0 [nm]
%
%     Seg1 rate = (w_L1_nm - w_start_nm) / L1_um  — fixed (L1 is constant).
%     Seg2 rate = (w_end_nm - w_L1_nm) / (L - L1_um)  — computed per-L at runtime
%     because seg2 length varies across the length sweep.
%
%     Cases without L1_um use a single dk profile (v3.2 behaviour, backward compat).
%     If L <= L1_um, only segment 1 is simulated.
%
%   - Phase and loss across the segment boundary:
%       dk_center_2 = dk_center_1 + dk_grad_1 * L1   (dk is continuous at boundary)
%       Phase: Phi(z_local) = Phi(L1) + dk_center_2*z_local + 0.5*dk_grad_2*z_local^2
%       State vector Y is continuous (no reset at boundary).
%       Loss gradient uses absolute z throughout so a(z) is smooth across segments.
%
% Inherited from v3.2 (all prior features unchanged):
%   - Z-dependent loss: a0_grad_dBcm_mm, a3_grad_dBcm_mm per case
%   - Single-segment geometry-based chirp via taper struct (dh_per_mm, dw_per_mm)
%   - Z-dependent phase mismatch: dk(z) = dk_center + dk_grad * z
%   - Accumulated phase Phi(z) = dk_center*z + 0.5*dk_grad*z^2
%   - Scattered SHG tracking via RK4 state Y(5)
%   - PM tuning curves, loss sweep, interactive evaluation
%
% RELATED SCRIPTS IN THIS FOLDER:
%   SHG_Design_Suite_v2_1.m    -- kept for reference only, not maintained. Last
%                                 version with a working numeric-vs-analytic
%                                 validation path (separate g_numeric/g_analytic,
%                                 run_validation_plot). Simpler: one case, no
%                                 z-dependent gradients/taper, no dispersion table.
%                                 Uses exp(+i*dk*z) phase convention.
%   SHG_Efficiency_w_Scatter.m -- a different kind of study: sweeps pump
%                                 wavelength detuning (not length/loss) to compare
%                                 guided vs. scattered SHG phase-matching bandwidth.
%                                 Shares this script's RK4 engine (helpers/shg/)
%                                 but keeps its own, less-precise g-coefficient
%                                 constants -- not directly comparable number-for-number.

    close all; clear; clc;
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\shg');

    %% --- 1. GLOBAL SIMULATION CONFIGURATION ---
    sim.dz = 0.25;   % RK4 step size (um)

    % --- SWEEP SETUP ---
    study_type  = 'length';
    sweep_vals  = linspace(100, 3000, 50);
    param_label = 'Waveguide Length [\mum]';

    % --- PHASE MATCHING TUNING ---
    dk_tune.L_vals_um  = [500, 1000, 2000, 3000];
    dk_tune.n_pts      = 201;
    dk_tune.n_lobes    = 3;
    dk_tune.dl_max_nm  = 5;
    dk_tune.dz_um      = 1.0;

    %% --- 2. CASE DEFINITIONS ---

    defaults.name             = '';
    defaults.Pp_avg_mW        = 0;
    defaults.duty_factor      = 1.0;
    defaults.lam_nm           = 458;
    defaults.d33_pmV          = 7.0;
    defaults.n_pump           = 2.05;
    defaults.n_shg            = 2.05;
    defaults.overlap_eta      = 0.0152;
    defaults.width_um         = 0.400;
    defaults.height_um        = 0.335;
    defaults.OCE              = 0.15;
    defaults.uv_loss_mode     = 'multiplier';
    defaults.uv_loss_val      = 10;
    defaults.dk_center        = 0;
    defaults.dk_grad          = 0;
    defaults.a0_grad_dBcm_mm  = 0;
    defaults.a3_grad_dBcm_mm  = 0;

    % --- CASE 1 ---
    cases(1) = defaults;
    cases(1).name        = 'Benchtop RunF (Avo#5)';
    cases(1).Pp_avg_mW   = 10 * 0.05 * 0.8;
    cases(1).duty_factor = 4e-5;

    % --- CASE 2 ---
    cases(2) = defaults;
    cases(2).name      = 'Gen 0 Pkg';
    cases(2).Pp_avg_mW = 200 * 0.6 * 0.8;

    % --- CASE 3 ---
    cases(3) = defaults;
    cases(3).name            = 'Pump Laser: 600mW';
    cases(3).Pp_avg_mW       = 600 * 0.6 * 0.8;
    cases(3).a0_grad_dBcm_mm = -10;

    % --- CASE 4 ---
    cases(4) = defaults;
    cases(4).name            = 'Input Coupling Improvement (60% -> 80%)';
    cases(4).Pp_avg_mW       = 600 * 0.8 * 0.8;
    cases(4).a0_grad_dBcm_mm = -10;

    % --- CASE 5 ---
    cases(5) = defaults;
    cases(5).name      = 'Input Facet ARC';
    cases(5).Pp_avg_mW = 600 * 0.8 * 1.0;

    % --- CASE 6 ---
    cases(6) = defaults;
    cases(6).name      = 'Pulsing (2x)';
    cases(6).Pp_avg_mW = 1200 * 0.8 * 1.0;

    % --- CASE 7: single-segment height taper (v3.2 reference) ---
    cases(7) = defaults;
    cases(7).name      = 'Output Coupler';
    cases(7).Pp_avg_mW = 1200 * 0.8 * 1.0;
    cases(7).OCE       = 0.85;

    % --- CASE 8: two-segment width taper ---
    cases(8) = defaults;
    cases(8).name      = 'Two-Segment Width Taper';
    cases(8).Pp_avg_mW = 1200 * 0.8 * 1.0;
    cases(8).OCE       = 0.85;
    % Specify widths [nm] at key points; L1_um activates segmented mode.
    cases(8).taper.lam_PM_0_nm  = 458.0;   % pump WL at which dk=0 for the cases(k).width_um geometry [nm]
    cases(8).taper.L1_um        = 350;      % segment 1 fixed length [um]
    cases(8).taper.w_start_nm   = 2500;      % width at z=0 (input) [nm]
    cases(8).taper.w_L1_nm      = 400;      % width at seg1/seg2 boundary [nm]
    cases(8).taper.w_end_nm     = 400;      % width at waveguide output [nm]
    cases(8).taper.dlam_PM_dw   = -0.012;  % width PM sensitivity [nm/nm]
    cases(8).taper.dh_per_mm    =  0;    % height taper: global [nm/mm]
    cases(8).taper.dlam_PM_dh   =  0.22;   % height PM sensitivity [nm/nm]

    %% --- 3. LOSS SWEEP CONFIGURATION ---
    blue_losses = [50, 35, 25, 20, 10, 5];
    num_losses  = length(blue_losses);
    num_cases   = length(cases);
    num_lengths = length(sweep_vals);

    dk_tune.loss_idx = [1, num_losses];

    %% --- 3b. RESOLVE DISPERSION AND TAPER PARAMETERS ---
    script_dir = fileparts(mfilename('fullpath'));
    json_path  = fullfile(script_dir, 'dispersion_table.json');
    [dk_tune.dn_pump_per_nm, dk_tune.dn_shg_per_nm] = ...
        shg_resolve_dispersion(cases(end), json_path);
    cases = shg_resolve_taper_params(cases, dk_tune);

    %% --- 4. EXECUTE SIMULATION ---
    master_guided = zeros(num_lengths, num_losses, num_cases);
    master_scat   = zeros(num_lengths, num_losses, num_cases);

    fprintf('SHG Design Suite v3.3\n');
    fprintf('=====================\n');

    for k = 1:num_cases
        fprintf('\nCase %d/%d: %s\n', k, num_cases, cases(k).name);

        cfg.sim  = sim;
        cfg.phys = cases(k);
        cfg.phys.Pp_peak_mW = cfg.phys.Pp_avg_mW / cfg.phys.duty_factor;
        cfg.phys.g          = shg_g_coefficient(cfg);

        A_eff_m2  = (cfg.phys.width_um * cfg.phys.height_um) * 1e-12;
        ol_factor = cfg.phys.overlap_eta / sqrt(A_eff_m2);
        fprintf('  overlap_eta = %.4f  (%.2f%%)    ol_factor = %.4ge6 /m    g = %.4e /sqrt(W)/um\n', ...
                cfg.phys.overlap_eta, cfg.phys.overlap_eta * 100, ol_factor / 1e6, cfg.phys.g);

        % Taper / dk summary
        if isfield(cfg.phys, 'L1_um') && ~isempty(cfg.phys.L1_um)
            dk_at_L1 = cfg.phys.dk_center + cfg.phys.dk_grad * cfg.phys.L1_um;
            fprintf('  Two-segment taper: L1 = %d um\n', cfg.phys.L1_um);
            fprintf('    Seg1: dk(0)=%.4e  dk(L1)=%.4e  dk_grad_1=%.3e rad/um^2\n', ...
                    cfg.phys.dk_center, dk_at_L1, cfg.phys.dk_grad);
            if isfield(cfg.phys, 'w_end_nm') && ~isempty(cfg.phys.w_end_nm)
                fprintf('    Seg2: w: %.0f nm -> %.0f nm  (dk_grad_2 computed per L)\n', ...
                        cfg.phys.w_L1_nm, cfg.phys.w_end_nm);
            else
                fprintf('    Seg2: dk_grad_2=%.3e rad/um^2\n', cfg.phys.dk_grad_2);
            end
        elseif cfg.phys.dk_grad ~= 0
            dk_end = cfg.phys.dk_center + cfg.phys.dk_grad * sweep_vals(end);
            fprintf('  dk_grad = %.3e rad/um^2   =>  dk(z=0) = %.4e,  dk(z=%.0fum) = %.4e rad/um\n', ...
                    cfg.phys.dk_grad, cfg.phys.dk_center, sweep_vals(end), dk_end);
        end

        if cfg.phys.a0_grad_dBcm_mm ~= 0 || cfg.phys.a3_grad_dBcm_mm ~= 0
            fprintf('  Loss gradients:  a0_grad = %.4f dB/cm/mm,  a3_grad = %.4f dB/cm/mm\n', ...
                    cfg.phys.a0_grad_dBcm_mm, cfg.phys.a3_grad_dBcm_mm);
        end

        loss_cases = struct('a0', {}, 'a3_scat', {}, 'a3_abs', {});
        for c = 1:num_losses
            loss_cases(c).a0 = blue_losses(c);
            if strcmp(cfg.phys.uv_loss_mode, 'multiplier')
                loss_cases(c).a3_scat = blue_losses(c) * cfg.phys.uv_loss_val;
            else
                loss_cases(c).a3_scat = cfg.phys.uv_loss_val;
            end
            loss_cases(c).a3_abs = 0;
        end

        for c = 1:num_losses
            case_cfg = cfg;
            case_cfg.phys.a0_dBcm      = loss_cases(c).a0;
            case_cfg.phys.a3_scat_dBcm = loss_cases(c).a3_scat;
            case_cfg.phys.a3_abs_dBcm  = loss_cases(c).a3_abs;

            if cfg.phys.a0_grad_dBcm_mm == 0 && cfg.phys.a3_grad_dBcm_mm == 0
                fprintf('    a0 = %d dB/cm,  a3 = %d dB/cm\n', ...
                        loss_cases(c).a0, loss_cases(c).a3_scat);
            else
                L_max_mm = sweep_vals(end) * 1e-3;
                fprintf('    a0 = %d dB/cm at z=0,  a0(z=%.0fum) = %.1f dB/cm  |  a3 = %d dB/cm at z=0\n', ...
                        loss_cases(c).a0, sweep_vals(end), ...
                        loss_cases(c).a0 + cfg.phys.a0_grad_dBcm_mm * L_max_mm, ...
                        loss_cases(c).a3_scat);
            end

            for i = 1:num_lengths
                if strcmp(study_type, 'length')
                    case_cfg.phys.L_um = sweep_vals(i);
                end
                [Pg_peak, Ps_peak] = shg_run_core_simulation(case_cfg);
                master_guided(i, c, k) = Pg_peak * cfg.phys.duty_factor * cfg.phys.OCE;
                master_scat(i, c, k)   = Ps_peak * cfg.phys.duty_factor * cfg.phys.OCE;
            end
        end

        shg_generate_case_plots(cases(k), sweep_vals, master_guided(:,:,k), master_scat(:,:,k), blue_losses, param_label);
    end

    %% --- 5. SHG vs LOSS PLOTS (fixed lengths) ---
    shg_generate_loss_axis_plot(cases, master_guided, master_scat, sweep_vals, blue_losses,  500);
    shg_generate_loss_axis_plot(cases, master_guided, master_scat, sweep_vals, blue_losses, 2000);

    %% --- 6. PHASE MATCHING TUNING CURVES ---
    fprintf('\nGenerating phase matching tuning curves...\n');
    shg_generate_pm_tuning_plot(cases, dk_tune, sim, blue_losses);

    %% --- 8. INTERACTIVE TARGET EVALUATION ---
    fprintf('\n===================================================\n');
    fprintf('           POST-RUN EVALUATION INTERFACE           \n');
    fprintf('===================================================\n');

    user_input = input('Enter evaluation target length in um [Default = 2000]: ', 's');
    if isempty(user_input)
        target_L = 2000;
    else
        target_L = str2double(user_input);
    end
    if isnan(target_L) || target_L < sweep_vals(1) || target_L > sweep_vals(end)
        fprintf('Invalid or out-of-bounds. Defaulting to 2000 um.\n');
        target_L = 2000;
    end

    matrix_guided = zeros(num_cases, num_losses);
    matrix_scat   = zeros(num_cases, num_losses);
    for k = 1:num_cases
        for c = 1:num_losses
            matrix_guided(k, c) = interp1(sweep_vals, master_guided(:, c, k), target_L, 'linear', 'extrap');
            matrix_scat(k, c)   = interp1(sweep_vals, master_scat(:, c, k),   target_L, 'linear', 'extrap');
        end
    end

    fprintf('\n>>> GUIDED SHG POWER  P_guided  (at L = %.0f um) <<<\n', target_L);
    shg_print_table(cases, blue_losses, matrix_guided, num_cases, num_losses);

    fprintf('\n>>> SCATTERED SHG POWER  P_scat  (at L = %.0f um) <<<\n', target_L);
    shg_print_table(cases, blue_losses, matrix_scat, num_cases, num_losses);

    fprintf('\n>>> TOTAL SHG POWER  P_guided + P_scat  (at L = %.0f um) <<<\n', target_L);
    shg_print_table(cases, blue_losses, matrix_guided + matrix_scat, num_cases, num_losses);

    csv_choice = input('Export guided SHG matrix to CSV? (y/n) [n]: ', 's');
    if strcmpi(csv_choice, 'y') || strcmpi(csv_choice, 'yes')
        filename = sprintf('SHG_v3_3_Evaluation_%.0fum.csv', target_L);
        fid = fopen(filename, 'w');
        fprintf(fid, 'Case ID,Case Name');
        for c = 1:num_losses
            fprintf(fid, ',Guided_Alpha0_%ddBcm_uW,Scat_Alpha0_%ddBcm_uW', blue_losses(c), blue_losses(c));
        end
        fprintf(fid, '\n');
        for k = 1:num_cases
            fprintf(fid, '%d,"%s"', k, cases(k).name);
            for c = 1:num_losses
                fprintf(fid, ',%.6f,%.6f', matrix_guided(k,c)*1e6, matrix_scat(k,c)*1e6);
            end
            fprintf(fid, '\n');
        end
        fclose(fid);
        fprintf('Saved: %s\n', filename);
    end
end

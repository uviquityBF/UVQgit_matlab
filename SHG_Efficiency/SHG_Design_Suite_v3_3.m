function SHG_Design_Suite_v3_3()
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

    close all; clear; clc;

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

    % --- CASE 1 ---
    cases(1).name             = 'Benchtop RunF (Avo#5)';
    cases(1).Pp_avg_mW        = 10 * 0.05 * 0.8;
    cases(1).duty_factor      = 4e-5;
    cases(1).lam_nm           = 450;
    cases(1).d33_pmV          = 7.0;
    cases(1).n_pump           = 2.05;
    cases(1).n_shg            = 2.05;
    cases(1).overlap_eta      = 0.0152;
    cases(1).width_um         = 0.300;
    cases(1).height_um        = 0.335;
    cases(1).OCE              = 0.15;
    cases(1).uv_loss_mode     = 'multiplier';
    cases(1).uv_loss_val      = 10;
    cases(1).dk_center        = 10;
    cases(1).dk_grad          = 0;
    cases(1).a0_grad_dBcm_mm  = 0;
    cases(1).a3_grad_dBcm_mm  = 0;

    % --- CASE 2 ---
    cases(2).name             = 'Gen 0 Pkg';
    cases(2).Pp_avg_mW        = 200 * 0.6 * 0.8;
    cases(2).duty_factor      = 1.0;
    cases(2).lam_nm           = 450;
    cases(2).d33_pmV          = 7.0;
    cases(2).n_pump           = 2.05;
    cases(2).n_shg            = 2.05;
    cases(2).overlap_eta      = 0.0152;
    cases(2).width_um         = 0.300;
    cases(2).height_um        = 0.335;
    cases(2).OCE              = 0.15;
    cases(2).uv_loss_mode     = 'multiplier';
    cases(2).uv_loss_val      = 10;
    cases(2).dk_center        = 0;
    cases(2).dk_grad          = 0;
    cases(2).a0_grad_dBcm_mm  = 0;
    cases(2).a3_grad_dBcm_mm  = 0;
    cases(2).taper.lam_PM_0_nm = 450.0;
    cases(2).taper.dlam_PM_dh  = 0.22;
    cases(2).taper.dlam_PM_dw  = -0.012;
    cases(2).taper.dh_per_mm   = 3.0;
    cases(2).taper.dw_per_mm   = 0.0;

    % --- CASE 3 ---
    cases(3).name             = 'Pump Laser: 600mW';
    cases(3).Pp_avg_mW        = 600 * 0.6 * 0.8;
    cases(3).duty_factor      = 1.0;
    cases(3).lam_nm           = 450;
    cases(3).d33_pmV          = 7.0;
    cases(3).n_pump           = 2.05;
    cases(3).n_shg            = 2.05;
    cases(3).overlap_eta      = 0.0152;
    cases(3).width_um         = 0.300;
    cases(3).height_um        = 0.335;
    cases(3).OCE              = 0.15;
    cases(3).uv_loss_mode     = 'multiplier';
    cases(3).uv_loss_val      = 10;
    cases(3).dk_center        = 0;
    cases(3).dk_grad          = 0;
    cases(3).a0_grad_dBcm_mm  = -10;
    cases(3).a3_grad_dBcm_mm  = 0;

    % --- CASE 4 ---
    cases(4).name             = 'Input Coupling Improvement (60% -> 80%)';
    cases(4).Pp_avg_mW        = 600 * 0.8 * 0.8;
    cases(4).duty_factor      = 1.0;
    cases(4).lam_nm           = 450;
    cases(4).d33_pmV          = 7.0;
    cases(4).n_pump           = 2.05;
    cases(4).n_shg            = 2.05;
    cases(4).overlap_eta      = 0.0152;
    cases(4).width_um         = 0.300;
    cases(4).height_um        = 0.335;
    cases(4).OCE              = 0.15;
    cases(4).uv_loss_mode     = 'multiplier';
    cases(4).uv_loss_val      = 10;
    cases(4).dk_center        = 0;
    cases(4).dk_grad          = 0;
    cases(4).a0_grad_dBcm_mm  = -10;
    cases(4).a3_grad_dBcm_mm  = 0;

    % --- CASE 5 ---
    cases(5).name             = 'Input Facet ARC';
    cases(5).Pp_avg_mW        = 600 * 0.8 * 1.0;
    cases(5).duty_factor      = 1.0;
    cases(5).lam_nm           = 450;
    cases(5).d33_pmV          = 7.0;
    cases(5).n_pump           = 2.05;
    cases(5).n_shg            = 2.05;
    cases(5).overlap_eta      = 0.0152;
    cases(5).width_um         = 0.300;
    cases(5).height_um        = 0.335;
    cases(5).OCE              = 0.15;
    cases(5).uv_loss_mode     = 'multiplier';
    cases(5).uv_loss_val      = 10;
    cases(5).dk_center        = 0;
    cases(5).dk_grad          = 0;
    cases(5).a0_grad_dBcm_mm  = 0;
    cases(5).a3_grad_dBcm_mm  = 0;

    % --- CASE 6 ---
    cases(6).name             = 'Pulsing (2x)';
    cases(6).Pp_avg_mW        = 1200 * 0.8 * 1.0;
    cases(6).duty_factor      = 1.0;
    cases(6).lam_nm           = 450;
    cases(6).d33_pmV          = 7.0;
    cases(6).n_pump           = 2.05;
    cases(6).n_shg            = 2.05;
    cases(6).overlap_eta      = 0.0152;
    cases(6).width_um         = 0.300;
    cases(6).height_um        = 0.335;
    cases(6).OCE              = 0.15;
    cases(6).uv_loss_mode     = 'multiplier';
    cases(6).uv_loss_val      = 10;
    cases(6).dk_center        = 0;
    cases(6).dk_grad          = 0;
    cases(6).a0_grad_dBcm_mm  = 0;
    cases(6).a3_grad_dBcm_mm  = 0;

    % --- CASE 7: single-segment height taper (v3.2 reference) ---
    cases(7).name             = 'Output Coupler';
    cases(7).Pp_avg_mW        = 1200 * 0.8 * 1.0;
    cases(7).duty_factor      = 1.0;
    cases(7).lam_nm           = 450;
    cases(7).d33_pmV          = 7.0;
    cases(7).n_pump           = 2.05;
    cases(7).n_shg            = 2.05;
    cases(7).overlap_eta      = 0.0152;
    cases(7).width_um         = 0.300;
    cases(7).height_um        = 0.335;
    cases(7).OCE              = 0.85;
    cases(7).uv_loss_mode     = 'multiplier';
    cases(7).uv_loss_val      = 10;
    cases(7).dk_center        = 0;
    cases(7).dk_grad          = 0;
    cases(7).a0_grad_dBcm_mm  = 0;
    cases(7).a3_grad_dBcm_mm  = 0;
    %cases(7).taper.lam_PM_0_nm = 450.0;
    %cases(7).taper.dlam_PM_dh  = 0.22;
    %cases(7).taper.dlam_PM_dw  = -0.012;
    %cases(7).taper.dh_per_mm   = 3.0;
    %cases(7).taper.dw_per_mm   = 0.0;

    % --- CASE 8: two-segment width taper (v3.3 new feature) ---
    % Segment 1 (L1=800um): combined height + width chirp
    % Segment 2 (remainder): height chirp only (width taper stops at L1)
    % This demonstrates the piecewise taper: seg1 varies width, seg2 continues
    % with only the global height taper.
    cases(8).name             = 'Two-Segment Width Taper';
    cases(8).Pp_avg_mW        = 1200 * 0.8 * 1.0;
    cases(8).duty_factor      = 1.0;
    cases(8).lam_nm           = 450;
    cases(8).d33_pmV          = 7.0;
    cases(8).n_pump           = 2.05;
    cases(8).n_shg            = 2.05;
    cases(8).overlap_eta      = 0.0152;
    cases(8).width_um         = 0.300;
    cases(8).height_um        = 0.335;
    cases(8).OCE              = 0.85;
    cases(8).uv_loss_mode     = 'multiplier';
    cases(8).uv_loss_val      = 10;
    cases(8).dk_center        = 0;
    cases(8).dk_grad          = 0;
    cases(8).a0_grad_dBcm_mm  = 0;
    cases(8).a3_grad_dBcm_mm  = 0;
    % Two-segment width taper: L1_um activates segmented mode.
    % Specify actual widths [nm] at three points — input, seg boundary, output.
    cases(8).taper.lam_PM_0_nm  = 450.0;   % PM wavelength at waveguide input [nm]
    cases(8).taper.L1_um        = 800;      % segment 1 fixed length [um]
    cases(8).taper.w_start_nm   = 300;      % width at z=0 (input) [nm]
    cases(8).taper.w_L1_nm      = 310;      % width at seg1/seg2 boundary [nm]
    cases(8).taper.w_end_nm     = 310;      % width at waveguide output [nm]
    cases(8).taper.dlam_PM_dw   = -0.012;  % width PM sensitivity [nm/nm]
    cases(8).taper.dh_per_mm    =  3.0;    % height taper: global [nm/mm]
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
        resolve_dispersion(cases(end), json_path);
    cases = resolve_taper_params(cases, dk_tune);

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
        cfg.phys.g          = calculate_g(cfg);

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
                [Pg_peak, Ps_peak] = run_core_simulation(case_cfg);
                master_guided(i, c, k) = Pg_peak * cfg.phys.duty_factor * cfg.phys.OCE;
                master_scat(i, c, k)   = Ps_peak * cfg.phys.duty_factor * cfg.phys.OCE;
            end
        end

        generate_case_plots(cases(k), sweep_vals, master_guided(:,:,k), master_scat(:,:,k), blue_losses, param_label);
    end

    %% --- 5. SHG vs LOSS PLOTS (fixed lengths) ---
    generate_loss_axis_plot(cases, master_guided, master_scat, sweep_vals, blue_losses,  500);
    generate_loss_axis_plot(cases, master_guided, master_scat, sweep_vals, blue_losses, 2000);

    %% --- 6. PHASE MATCHING TUNING CURVES ---
    fprintf('\nGenerating phase matching tuning curves...\n');
    generate_pm_tuning_plot(cases, dk_tune, sim, blue_losses);

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
    print_table(cases, blue_losses, matrix_guided, num_cases, num_losses);

    fprintf('\n>>> SCATTERED SHG POWER  P_scat  (at L = %.0f um) <<<\n', target_L);
    print_table(cases, blue_losses, matrix_scat, num_cases, num_losses);

    fprintf('\n>>> TOTAL SHG POWER  P_guided + P_scat  (at L = %.0f um) <<<\n', target_L);
    print_table(cases, blue_losses, matrix_guided + matrix_scat, num_cases, num_losses);

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


%% =========================================================================
%  G-COEFFICIENT  (Wei's formula)
% =========================================================================
function g = calculate_g(cfg)
    eps0   = 8.854e-12;
    c      = 2.998e8;
    lam1   = cfg.phys.lam_nm * 1e-9;
    omega1 = 2 * pi * c / lam1;
    deff   = cfg.phys.d33_pmV * 1e-12;
    n1     = cfg.phys.n_pump;
    n2     = cfg.phys.n_shg;
    A_eff  = (cfg.phys.width_um * cfg.phys.height_um) * 1e-12;

    ol_factor = cfg.phys.overlap_eta / sqrt(A_eff);
    yeta      = 2 * omega1^2 / (eps0 * c^3 * n1^2 * n2) * deff^2 * ol_factor^2;
    g         = sqrt(yeta) * 1e-6;
end


%% =========================================================================
%  CORE SOLVER WRAPPER
%  Dispatches to single-segment or two-segment RK4 engine based on whether
%  L1_um is present and non-empty (set by resolve_taper_params for cases
%  that have a two-segment taper struct).
% =========================================================================
function [Pg, Ps] = run_core_simulation(cfg)
    Pp  = cfg.phys.Pp_peak_mW * 1e-3;
    a0  = (cfg.phys.a0_dBcm      / 4.3429) * 1e-4;
    as  = (cfg.phys.a3_scat_dBcm / 4.3429) * 1e-4;
    aa  = (cfg.phys.a3_abs_dBcm  / 4.3429) * 1e-4;
    a3t = as + aa;

    a0_grad  = (cfg.phys.a0_grad_dBcm_mm / 4.3429) * 1e-7;
    a3t_grad = (cfg.phys.a3_grad_dBcm_mm / 4.3429) * 1e-7;
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

        [Pg, Ps] = rk4_engine_segmented( ...
            cfg.phys.L_um, L1, cfg.sim.dz, cfg.phys.g, ...
            a0, a0_grad, a3t, a3t_grad, as, a3s_grad, ...
            cfg.phys.dk_center, cfg.phys.dk_grad, dk_grad_2, Pp);
    else
        [Pg, Ps] = rk4_engine( ...
            cfg.phys.L_um, cfg.sim.dz, cfg.phys.g, ...
            a0, a0_grad, a3t, a3t_grad, as, a3s_grad, ...
            cfg.phys.dk_center, cfg.phys.dk_grad, Pp);
    end
end


%% =========================================================================
%  RK4 ENGINE — single segment (backward-compatible wrapper)
% =========================================================================
function [Pg_end, Ps_end] = rk4_engine(L, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, dk_center, dk_grad, Pp)
    Y0 = [sqrt(Pp); 0; 0; 0; 0];
    [Y_end, ~] = rk4_segment(L, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, ...
                              dk_center, dk_grad, Y0, 0.0, 0.0);
    Pg_end = Y_end(3)^2 + Y_end(4)^2;
    Ps_end = Y_end(5);
end


%% =========================================================================
%  RK4 ENGINE — two-segment
%
%  seg1: z_local 0 -> L1,   z_abs 0 -> L1,   phase starts at 0
%  seg2: z_local 0 -> L-L1, z_abs L1 -> L,   phase starts at phase_end_seg1
%
%  dk_center_2 = dk_center_1 + dk_grad_1 * L1  (dk is continuous at boundary)
%
%  If L <= L1, only seg1 runs (waveguide shorter than the first segment).
% =========================================================================
function [Pg_end, Ps_end] = rk4_engine_segmented(L, L1, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, dk_center_1, dk_grad_1, dk_grad_2, Pp)
    Y0 = [sqrt(Pp); 0; 0; 0; 0];

    if L <= L1
        [Y_end, ~] = rk4_segment(L, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, ...
                                  dk_center_1, dk_grad_1, Y0, 0.0, 0.0);
    else
        [Y1, phase1] = rk4_segment(L1, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, ...
                                    dk_center_1, dk_grad_1, Y0, 0.0, 0.0);
        dk_center_2 = dk_center_1 + dk_grad_1 * L1;
        L2 = L - L1;
        [Y_end, ~] = rk4_segment(L2, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, ...
                                  dk_center_2, dk_grad_2, Y1, L1, phase1);
    end

    Pg_end = Y_end(3)^2 + Y_end(4)^2;
    Ps_end = Y_end(5);
end


%% =========================================================================
%  RK4 SEGMENT — core loop for one piecewise section
%
%  z_loc  : position local to this segment (0 at segment start)  — drives phase
%  z_abs  : absolute position from waveguide input = z_loc + z_offset — drives loss
%  phase  : total accumulated phase = phase_offset + dk_center*z_loc + 0.5*dk_grad*z_loc^2
%
%  Returns Y_end (state at end of segment) and phase_end (total accumulated
%  phase at end of segment, to be passed as phase_offset to the next segment).
% =========================================================================
function [Y_end, phase_end] = rk4_segment(L_seg, dz, g, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, dk_center, dk_grad, Y0, z_offset, phase_offset)
    z_vec = 0:dz:L_seg;
    Y     = Y0;

    for i = 1:(length(z_vec) - 1)
        z_loc = z_vec(i);
        z_abs = z_loc + z_offset;
        k1 = shg_derivs_seg(z_loc,       z_abs,       Y,           a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset);
        k2 = shg_derivs_seg(z_loc+dz/2, z_abs+dz/2, Y+k1*dz/2,   a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset);
        k3 = shg_derivs_seg(z_loc+dz/2, z_abs+dz/2, Y+k2*dz/2,   a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset);
        k4 = shg_derivs_seg(z_loc+dz,   z_abs+dz,   Y+k3*dz,     a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset);
        Y  = Y + (dz/6) * (k1 + 2*k2 + 2*k3 + k4);
    end

    Y_end     = Y;
    phase_end = phase_offset + dk_center * L_seg + 0.5 * dk_grad * L_seg^2;
end


%% =========================================================================
%  CME DERIVATIVES — one segment
%
%  z_loc        : local z within this segment [um]  — used for phase integral
%  z_abs        : absolute z from waveguide input [um]  — used for loss gradient
%  phase_offset : phase accumulated in all prior segments [rad]
%
%  Phase: Phi = phase_offset + dk_center*z_loc + 0.5*dk_grad*z_loc^2
%  Loss:  a(z) = a_center + a_grad * z_abs   (continuous across segments)
% =========================================================================
function dY = shg_derivs_seg(z_loc, z_abs, Y, a0, a0_grad, a3t, a3t_grad, a3s, a3s_grad, g, dk_center, dk_grad, phase_offset)
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


%% =========================================================================
%  RESOLVE TAPER PARAMETERS
%
%  Single-segment (no L1_um in taper struct):
%    dlam_PM_dz = (dlam_PM_dh*dh_per_mm + dlam_PM_dw*dw_per_mm) * 1e-3  [nm/um]
%    sets cases(k).dk_center and cases(k).dk_grad
%
%  Two-segment, width-endpoint spec (preferred): L1_um + w_start_nm/w_L1_nm/w_end_nm
%    Seg1 rate derived from (w_L1_nm - w_start_nm) / L1_um — fixed since L1 is constant.
%    Seg2 rate deferred: stored as dk_grad_h + C_scale_dw factors so run_core_simulation
%    can recompute dk_grad_2 per-L (seg2 length = L - L1 varies across the sweep).
%    Sets: cases(k).dk_center, dk_grad (seg1), L1_um, w_L1_nm, w_end_nm, dk_grad_h, C_scale_dw.
%
%  Two-segment, rate spec (legacy): L1_um + dw_per_mm_1 + dw_per_mm_2
%    Both dk_grad values pre-computed; sets cases(k).dk_grad_2 directly.
% =========================================================================
function cases = resolve_taper_params(cases, dk_tune)
    any_taper = false;
    for k = 1:length(cases)
        if ~isfield(cases(k), 'taper') || isempty(cases(k).taper)
            continue;
        end
        t = cases(k).taper;
        if ~isfield(t, 'lam_PM_0_nm') || isnan(t.lam_PM_0_nm)
            continue;
        end
        if ~any_taper
            fprintf('\nResolving geometry-based chirp parameters...\n');
            any_taper = true;
        end

        lam_nm    = cases(k).lam_nm;
        C_rad_nm2 = (4*pi/lam_nm) * (dk_tune.dn_shg_per_nm/2 - dk_tune.dn_pump_per_nm);
        C_scale   = C_rad_nm2 * 1000;   % rad/(um*nm)

        dk_ctr = C_scale * (lam_nm - t.lam_PM_0_nm);

        if isfield(t, 'L1_um') && ~isempty(t.L1_um)
            % --- Two-segment width taper ---
            dlam_dz_h = t.dlam_PM_dh * t.dh_per_mm * 1e-3;   % height contribution [nm/um], global

            if isfield(t, 'w_start_nm')
                % Width-endpoint spec: actual widths at input, boundary, and output.
                % Seg1 rate is fixed (L1 is constant); seg2 rate is deferred to run_core_simulation.
                dw_mm_1   = (t.w_L1_nm - t.w_start_nm) / (t.L1_um * 1e-3);   % [nm/mm]
                dlam_dz_1 = dlam_dz_h + t.dlam_PM_dw * dw_mm_1 * 1e-3;
                dk_grd_1  = -C_scale * dlam_dz_1;

                cases(k).dk_center  = dk_ctr;
                cases(k).dk_grad    = dk_grd_1;
                cases(k).L1_um      = t.L1_um;
                cases(k).w_L1_nm    = t.w_L1_nm;
                cases(k).w_end_nm   = t.w_end_nm;
                cases(k).dk_grad_h  = -C_scale * dlam_dz_h;         % height-only dk gradient [rad/um^2]
                cases(k).C_scale_dw = -C_scale * t.dlam_PM_dw * 1e-3;  % width factor [rad/um^2 per nm/mm]

                fprintf('  Case %d (%s):  [two-segment, width endpoint]\n', k, cases(k).name);
                fprintf('    lam_PM_0 = %.2f nm  =>  dk_center = %.4e rad/um\n', t.lam_PM_0_nm, dk_ctr);
                fprintf('    Seg1 (L1=%d um): w %.0f->%.0f nm  dw=%.2f nm/mm  =>  dk_grad_1 = %.4e rad/um^2\n', ...
                        t.L1_um, t.w_start_nm, t.w_L1_nm, dw_mm_1, dk_grd_1);
                fprintf('    Seg2: w %.0f->%.0f nm  (dk_grad_2 computed per L at runtime)\n', ...
                        t.w_L1_nm, t.w_end_nm);
                fprintf('    C_scale = %.4e rad/(um*nm)\n', C_scale);
            else
                % Legacy rate spec: dw_per_mm_1 and dw_per_mm_2 (pre-computable).
                dlam_dz_1 = dlam_dz_h + t.dlam_PM_dw * t.dw_per_mm_1 * 1e-3;
                dlam_dz_2 = dlam_dz_h + t.dlam_PM_dw * t.dw_per_mm_2 * 1e-3;
                dk_grd_1  = -C_scale * dlam_dz_1;
                dk_grd_2  = -C_scale * dlam_dz_2;

                cases(k).dk_center = dk_ctr;
                cases(k).dk_grad   = dk_grd_1;
                cases(k).dk_grad_2 = dk_grd_2;
                cases(k).L1_um     = t.L1_um;

                fprintf('  Case %d (%s):  [two-segment, rate spec]\n', k, cases(k).name);
                fprintf('    lam_PM_0 = %.2f nm  =>  dk_center = %.4e rad/um\n', t.lam_PM_0_nm, dk_ctr);
                fprintf('    Seg1 (L1=%d um):   dlam_PM/dz = %.4f nm/mm  =>  dk_grad_1 = %.4e rad/um^2\n', t.L1_um, dlam_dz_1*1e3, dk_grd_1);
                fprintf('    Seg2 (remainder): dlam_PM/dz = %.4f nm/mm  =>  dk_grad_2 = %.4e rad/um^2\n', dlam_dz_2*1e3, dk_grd_2);
                fprintf('    C_scale = %.4e rad/(um*nm)\n', C_scale);
            end
        else
            % --- Single segment (v3.2 behavior) ---
            dlam_PM_dz = (t.dlam_PM_dh * t.dh_per_mm + t.dlam_PM_dw * t.dw_per_mm) * 1e-3;
            dk_grd     = -C_scale * dlam_PM_dz;

            cases(k).dk_center = dk_ctr;
            cases(k).dk_grad   = dk_grd;

            fprintf('  Case %d (%s):\n', k, cases(k).name);
            fprintf('    lam_PM_0 = %.2f nm  =>  dk_center = %.4e rad/um\n', t.lam_PM_0_nm, dk_ctr);
            fprintf('    dlam_PM/dz = %.4f nm/mm  =>  dk_grad = %.4e rad/um^2\n', dlam_PM_dz*1e3, dk_grd);
            fprintf('    C_scale = %.4e rad/(um*nm)\n', C_scale);
        end
    end
end


%% =========================================================================
%  PHASE MATCHING TUNING CURVE
% =========================================================================
function generate_pm_tuning_plot(cases, dk_tune, sim, blue_losses)
    ref_case = cases(end);

    lam_nm    = ref_case.lam_nm;
    dn_pump   = dk_tune.dn_pump_per_nm;
    dn_shg    = dk_tune.dn_shg_per_nm;
    C_rad_nm2 = (4*pi/lam_nm) * (dn_shg/2 - dn_pump);
    C_scale   = C_rad_nm2 * 1000;

    [a0_sorted, ~] = sort(blue_losses(dk_tune.loss_idx), 'ascend');
    n_loss_sel     = length(a0_sorted);

    cfg.sim              = sim;
    cfg.sim.dz           = dk_tune.dz_um;
    cfg.phys             = ref_case;
    cfg.phys.Pp_peak_mW  = cfg.phys.Pp_avg_mW / cfg.phys.duty_factor;
    cfg.phys.g           = calculate_g(cfg);
    cfg.phys.a3_abs_dBcm = 0;

    L_vals  = dk_tune.L_vals_um;
    n_L     = length(L_vals);
    n_pts   = dk_tune.n_pts;
    colors  = lines(n_L);
    dk_max  = abs(C_scale) * dk_tune.dl_max_nm;
    dk_vec  = linspace(-dk_max, dk_max, n_pts);

    all_dl_nm  = cell(n_L, n_loss_sel);
    all_guided = cell(n_L, n_loss_sel);
    all_scat   = cell(n_L, n_loss_sel);
    fwhm_g_nm  = NaN(n_L, n_loss_sel);
    fwhm_s_nm  = NaN(n_L, n_loss_sel);

    for lj = 1:n_loss_sel
        a0_ref = a0_sorted(lj);
        cfg.phys.a0_dBcm      = a0_ref;
        cfg.phys.a3_scat_dBcm = a0_ref * ref_case.uv_loss_val;
        fprintf('\n  Loss level: a0=%d dB/cm\n', a0_ref);

        for li = 1:n_L
            L = L_vals(li);
            guided_raw = zeros(1, n_pts);
            scat_raw   = zeros(1, n_pts);
            fprintf('    L=%d um (%d points)...\n', L, n_pts);

            for di = 1:n_pts
                cfg.phys.L_um      = L;
                cfg.phys.dk_center = dk_vec(di);
                [Pg, Ps] = run_core_simulation(cfg);
                guided_raw(di) = Pg * ref_case.duty_factor * ref_case.OCE;
                scat_raw(di)   = Ps * ref_case.duty_factor * ref_case.OCE;
            end

            pk_g  = max(guided_raw);
            pk_s  = max(scat_raw);
            dl_nm = dk_vec / C_scale;

            all_dl_nm{li, lj}  = dl_nm;
            all_guided{li, lj} = guided_raw / pk_g;
            all_scat{li, lj}   = scat_raw   / pk_s;
            fwhm_g_nm(li, lj)  = compute_pm_fwhm(dl_nm, all_guided{li, lj});
            fwhm_s_nm(li, lj)  = compute_pm_fwhm(dl_nm, all_scat{li, lj});
        end
    end

    disp_note = sprintf('dn_{pump}/d\\lambda=%.5f /nm,  dn_{SH}/d\\lambda=%.5f /nm', dn_pump, dn_shg);

    figure('Color','w','Name','Phase Matching Tuning Curves', ...
           'Position',[250, 60, 920, 320*n_loss_sel]);

    for lj = 1:n_loss_sel
        ax = subplot(n_loss_sel, 1, lj);
        hold(ax, 'on');

        h_leg      = gobjects(n_L, 1);
        leg_labels = cell(n_L, 1);

        for li = 1:n_L
            h_leg(li) = plot(ax, all_dl_nm{li,lj}, all_guided{li,lj}, '-', ...
                             'Color', colors(li,:), 'LineWidth', 1.8);
                         plot(ax, all_dl_nm{li,lj}, all_scat{li,lj},   '--', ...
                             'Color', colors(li,:), 'LineWidth', 1.2);

            g_fwhm = fwhm_g_nm(li, lj) * 1e3;
            s_fwhm = fwhm_s_nm(li, lj) * 1e3;
            if ~isnan(g_fwhm)
                leg_labels{li} = sprintf('L=%d \\mum  (G:%.0fpm  S:%.0fpm)', L_vals(li), g_fwhm, s_fwhm);
            else
                leg_labels{li} = sprintf('L=%d \\mum', L_vals(li));
            end
        end

        set(ax, 'YLim', [0 1.1]);  grid(ax, 'on');
        ylabel(ax, 'Norm. SHG Efficiency', 'FontSize', 10);
        title(ax, sprintf('\\alpha_0 = %d dB/cm  (SH = %d dB/cm)  —  %s  —  %s', ...
                           a0_sorted(lj), a0_sorted(lj)*ref_case.uv_loss_val, ...
                           ref_case.name, disp_note), 'FontSize', 8);
        legend(ax, h_leg, leg_labels, 'Location', 'northeast', 'FontSize', 8);
        text(0.02, 0.08, 'Solid = Guided,  Dashed = Scattered', ...
             'Units','normalized','FontSize',8,'Color',[0.35 0.35 0.35],'Parent',ax);

        if lj == n_loss_sel
            xlabel(ax, '\delta\lambda_{pump} (nm)  [detuning from perfect phase match;  \delta\lambda_{SHG} = \delta\lambda_{pump}/2]', 'FontSize', 10);
        end
    end

    fprintf('\nPhase Matching FWHM  (ref: %s)\n', ref_case.name);
    hdr = sprintf('%-10s', 'L (um)');
    for lj = 1:n_loss_sel
        hdr = [hdr sprintf('  | a0=%2ddB  G-FWHM(pm)  S-FWHM(pm)  S/G', a0_sorted(lj))]; %#ok
    end
    fprintf('%s\n%s\n', hdr, repmat('-', length(hdr)+2, 1));
    for li = 1:n_L
        row = sprintf('%-10d', L_vals(li));
        for lj = 1:n_loss_sel
            ratio = fwhm_s_nm(li,lj) / fwhm_g_nm(li,lj);
            row = [row sprintf('  |        %6.0f       %6.0f  %4.2f', ...
                               fwhm_g_nm(li,lj)*1e3, fwhm_s_nm(li,lj)*1e3, ratio)]; %#ok
        end
        fprintf('%s\n', row);
    end
    fprintf('%s\n', repmat('-', length(hdr)+2, 1));
end


%% =========================================================================
%  RESOLVE DISPERSION
% =========================================================================
function [dn_pump, dn_shg] = resolve_dispersion(ref_case, json_path)
    MEDIAN_PUMP = -0.001476;
    MEDIAN_SHG  = -0.013490;

    h_nm = ref_case.height_um * 1e3;
    w_nm = ref_case.width_um  * 1e3;

    if exist(json_path, 'file')
        try
            raw  = jsondecode(fileread(json_path));
            data = raw.data;
            h_t  = double([data.h_core_nm]');
            w_t  = double([data.w_core_nm]');
            dp_t = double([data.dn_dWL_fund]');
            ds_t = double([data.dn_dWL_target]');

            Fp = scatteredInterpolant(h_t, w_t, dp_t, 'linear', 'linear');
            Fs = scatteredInterpolant(h_t, w_t, ds_t, 'linear', 'linear');
            dn_pump = Fp(h_nm, w_nm);
            dn_shg  = Fs(h_nm, w_nm);

            in_bounds = check_table_bounds(h_t, w_t, h_nm, w_nm);
            fprintf('\nDispersion resolved from JSON  (h=%.0fnm, w=%.0fnm):\n', h_nm, w_nm);
            fprintf('  dn_pump/d\x03BB = %.6f /nm\n', dn_pump);
            fprintf('  dn_SH/d\x03BB   = %.6f /nm\n', dn_shg);
            if ~in_bounds
                fprintf('  WARNING: geometry outside table range (h=335-365nm, w=350-550nm) — extrapolated.\n');
            end
            return;
        catch err
            fprintf('\nWARNING: dispersion_table.json found but could not be read: %s\n', err.message);
        end
    else
        fprintf('\ndispersion_table.json not found in script directory.\n');
    end

    fprintf('Hardcoded median values (COMSOL table, 21 geometries):\n');
    fprintf('  dn_pump/d\x03BB = %.6f /nm\n', MEDIAN_PUMP);
    fprintf('  dn_SH/d\x03BB   = %.6f /nm\n', MEDIAN_SHG);

    resp = input('Accept these values? [Y/n]: ', 's');
    if strcmpi(strtrim(resp), 'n') || strcmpi(strtrim(resp), 'no')
        val = input(sprintf('  Enter dn_pump/d\x03BB [1/nm, default %.6f]: ', MEDIAN_PUMP), 's');
        dn_pump = str2double(strtrim(val));
        if isnan(dn_pump), dn_pump = MEDIAN_PUMP; end

        val = input(sprintf('  Enter dn_SH/d\x03BB   [1/nm, default %.6f]: ', MEDIAN_SHG), 's');
        dn_shg = str2double(strtrim(val));
        if isnan(dn_shg), dn_shg = MEDIAN_SHG; end
    else
        dn_pump = MEDIAN_PUMP;
        dn_shg  = MEDIAN_SHG;
    end

    fprintf('Using: dn_pump/d\x03BB = %.6f /nm,  dn_SH/d\x03BB = %.6f /nm\n', dn_pump, dn_shg);
end


%% =========================================================================
%  CHECK TABLE BOUNDS
% =========================================================================
function in_bounds = check_table_bounds(h_t, w_t, h, w)
    in_bounds = (h >= min(h_t)) && (h <= max(h_t)) && ...
                (w >= min(w_t)) && (w <= max(w_t));
end


%% =========================================================================
%  PM FWHM
% =========================================================================
function fwhm_dk = compute_pm_fwhm(dk_vec, shg_norm)
    fwhm_dk = NaN;
    [~, peak_idx] = max(shg_norm);
    half_max = 0.5;

    t_left = NaN;
    for i = peak_idx:-1:2
        if shg_norm(i-1) <= half_max && shg_norm(i) >= half_max
            t_left = interp1([shg_norm(i-1), shg_norm(i)], [dk_vec(i-1), dk_vec(i)], half_max);
            break;
        end
    end

    t_right = NaN;
    for i = peak_idx:(length(shg_norm)-1)
        if shg_norm(i) >= half_max && shg_norm(i+1) <= half_max
            t_right = interp1([shg_norm(i), shg_norm(i+1)], [dk_vec(i), dk_vec(i+1)], half_max);
            break;
        end
    end

    if ~isnan(t_left) && ~isnan(t_right)
        fwhm_dk = t_right - t_left;
    end
end


%% =========================================================================
%  PER-CASE PLOT
% =========================================================================
function generate_case_plots(case_struct, sweep_vals, guided, scat, blue_losses, param_label)
    num_losses    = length(blue_losses);
    colors        = lines(num_losses);
    h_guided      = gobjects(num_losses, 1);
    legend_labels = cell(num_losses, 1);

    for c = 1:num_losses
        legend_labels{c} = sprintf('\\alpha_0=%d dB/cm', blue_losses(c));
    end

    fig_title = sprintf('Case: %s', case_struct.name);
    sub_text  = sprintf('P_{avg}=%.1fmW, DF=%.2g (P_{peak}=%.1fmW), UV loss: %s x%g', ...
                case_struct.Pp_avg_mW, case_struct.duty_factor, ...
                case_struct.Pp_avg_mW / case_struct.duty_factor, ...
                case_struct.uv_loss_mode, case_struct.uv_loss_val);

    figure('Color','w','Name',fig_title,'Position',[150, 150, 870, 680]);

    subplot(2,1,1); hold on;
    for c = 1:num_losses
        h_guided(c) = semilogy(sweep_vals, guided(:,c)*1e6, '-', 'Color', colors(c,:), 'LineWidth', 1.8);
                      semilogy(sweep_vals, scat(:,c)*1e6,   '--','Color', colors(c,:), 'LineWidth', 1.2);
    end
    set(gca, 'YScale', 'log'); grid on;
    ylabel('Time-Avg Power [\muW]');
    title({fig_title; sub_text}, 'FontSize', 9);
    legend(h_guided, legend_labels, 'Location', 'eastoutside', 'FontSize', 8);
    text(0.02, 0.92, 'Solid = Guided,  Dashed = Scattered', ...
         'Units','normalized','FontSize',8,'Color',[0.3 0.3 0.3]);

    subplot(2,1,2); hold on;
    Pp_avg_W = case_struct.Pp_avg_mW * 1e-3;
    h2 = gobjects(num_losses, 1);
    for c = 1:num_losses
        Eff_g = (guided(:,c) ./ Pp_avg_W^2) * 100;
        Eff_s = (scat(:,c)   ./ Pp_avg_W^2) * 100;
        h2(c) = semilogy(sweep_vals, Eff_g, '-', 'Color', colors(c,:), 'LineWidth', 1.8);
                semilogy(sweep_vals, Eff_s, '--','Color', colors(c,:), 'LineWidth', 1.2);
    end
    set(gca, 'YScale', 'log'); grid on;
    ylabel('Normalized Efficiency [%/W^2]'); xlabel(param_label);
    title('Time-Averaged Intrinsic Conversion Efficiency', 'FontSize', 9);
    legend(h2, legend_labels, 'Location', 'eastoutside', 'FontSize', 8);
end


%% =========================================================================
%  TABLE PRINTER
% =========================================================================
function print_table(cases, blue_losses, matrix_W, num_cases, num_losses)
    fprintf('%-36s', 'Configuration');
    for c = 1:num_losses
        fprintf(' | a0=%2ddB/cm', blue_losses(c));
    end
    fprintf('\n%s\n', repmat('-', 38 + num_losses * 13, 1));
    for k = 1:num_cases
        fprintf('%-36s', sprintf('%d: %s', k, cases(k).name));
        for c = 1:num_losses
            fprintf(' | %10.3f uW', matrix_W(k, c) * 1e6);
        end
        fprintf('\n');
    end
    fprintf('%s\n', repmat('-', 38 + num_losses * 13, 1));
end


%% =========================================================================
%  SHG vs LOSS PLOT  (fixed waveguide length)
% =========================================================================
function generate_loss_axis_plot(cases, master_guided, master_scat, sweep_vals, blue_losses, target_L_um)
    num_cases  = length(cases);
    num_losses = length(blue_losses);
    colors     = lines(num_cases);

    guided_at_L = zeros(num_cases, num_losses);
    scat_at_L   = zeros(num_cases, num_losses);
    for k = 1:num_cases
        for c = 1:num_losses
            guided_at_L(k, c) = interp1(sweep_vals, master_guided(:, c, k), target_L_um, 'linear', 'extrap');
            scat_at_L(k, c)   = interp1(sweep_vals, master_scat(:, c, k),   target_L_um, 'linear', 'extrap');
        end
    end

    [loss_sorted, sort_idx] = sort(blue_losses);
    guided_sorted = guided_at_L(:, sort_idx);
    scat_sorted   = scat_at_L(:,   sort_idx);

    fig_name = sprintf('SHG vs Loss at L = %d um', target_L_um);
    figure('Color','w','Name', fig_name,'Position',[200, 200, 820, 560]);
    hold on;

    h = gobjects(num_cases, 1);
    for k = 1:num_cases
        h(k) = semilogy(loss_sorted, guided_sorted(k,:)*1e6, '-o', ...
                        'Color', colors(k,:), 'LineWidth', 1.8, 'MarkerSize', 6, ...
                        'MarkerFaceColor', colors(k,:));
               semilogy(loss_sorted, scat_sorted(k,:)*1e6, '--o', ...
                        'Color', colors(k,:), 'LineWidth', 1.2, 'MarkerSize', 6, ...
                        'MarkerFaceColor', 'w');
    end

    set(gca, 'YScale', 'log'); grid on;
    xlabel('Pump Loss \alpha_0 (dB/cm)', 'FontSize', 11);
    ylabel('Time-Avg SHG Power (\muW)', 'FontSize', 11);
    title({sprintf('SHG Power vs Waveguide Loss  —  L = %d \\mum', target_L_um); ...
           'Solid filled = Guided,  Dashed open = Scattered'}, 'FontSize', 10);

    case_labels = cell(num_cases, 1);
    for k = 1:num_cases
        case_labels{k} = sprintf('%d: %s', k, cases(k).name);
    end
    legend(h, case_labels, 'Location', 'eastoutside', 'FontSize', 7);
end

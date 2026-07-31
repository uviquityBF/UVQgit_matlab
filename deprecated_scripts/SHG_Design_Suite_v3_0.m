function SHG_Design_Suite_v3_0()
% SHG Design Suite v3.0
%
% Changes from v2_2:
%   - Refractive indices: n_pump = n_shg = 2.05 (Wei's effective index for both modes)
%   - Overlap input: overlap_eta (dimensionless, user-facing) converted internally to
%     ol_factor (1/m) via  ol_factor = overlap_eta / sqrt(A_eff).
%     Default value 0.0152 is Wei's ol_factor=0.048e6/m back-converted to dimensionless.
%   - g uses Wei's formula: g = sqrt( 2*w1^2*deff^2*ol_factor^2 / (eps0*c^3*n1^2*n2) )
%     This is algebraically identical to v2.2's g_numeric when ol_factor = eta/sqrt(A_eff).
%   - Phase sign in SH equation: exp(-j*dk*z)  (Wei's convention; no effect at dk=0)
%   - g_analytic removed (was computed but never used)
%   - g and overlap values printed at runtime for verification against Wei

    close all; clear; clc;

    %% --- 1. GLOBAL SIMULATION CONFIGURATION ---
    sim.dz = 0.25;   % RK4 step size (um)

    % --- SWEEP SETUP ---
    study_type  = 'length';
    sweep_vals  = linspace(100, 3000, 50);
    param_label = 'Waveguide Length [\mum]';

    %% --- 2. CASE DEFINITIONS ---
    %
    % overlap_eta: dimensionless mode overlap integral (COMSOL-style normalization).
    %   Internally converted to ol_factor (1/m) = overlap_eta / sqrt(width_um*height_um*1e-12).
    %   Wei reference: ol_factor = 0.048e6 /m  =>  overlap_eta ~ 0.0152  (~1.52%)
    %   for W=300nm, H=335nm waveguide, TM00+TM04 modes.
    %
    % n_pump = n_shg = 2.05  (Wei's effective index — same for both modes at this geometry)

    % --- CASE 1 ---
    cases(1).name           = 'Benchtop RunF (Avo#5)';
    cases(1).Pp_avg_mW      = 10 * 0.05 * 0.8;
    cases(1).duty_factor    = 4e-5;
    cases(1).lam_nm         = 450;
    cases(1).d33_pmV        = 7.0;
    cases(1).n_pump         = 2.05;
    cases(1).n_shg          = 2.05;
    cases(1).overlap_eta    = 0.0152;
    cases(1).width_um       = 0.300;
    cases(1).height_um      = 0.335;
    cases(1).OCE            = 0.15;
    cases(1).uv_loss_mode   = 'multiplier';
    cases(1).uv_loss_val    = 10;
    cases(1).dk_center      = 0;

    % --- CASE 2 ---
    cases(2).name           = 'Gen 0 Pkg';
    cases(2).Pp_avg_mW      = 200 * 0.6 * 0.8;
    cases(2).duty_factor    = 1.0;
    cases(2).lam_nm         = 450;
    cases(2).d33_pmV        = 7.0;
    cases(2).n_pump         = 2.05;
    cases(2).n_shg          = 2.05;
    cases(2).overlap_eta    = 0.0152;
    cases(2).width_um       = 0.300;
    cases(2).height_um      = 0.335;
    cases(2).OCE            = 0.15;
    cases(2).uv_loss_mode   = 'multiplier';
    cases(2).uv_loss_val    = 10;
    cases(2).dk_center      = 0;

    % --- CASE 3 ---
    cases(3).name           = 'Pump Laser: 600mW';
    cases(3).Pp_avg_mW      = 600 * 0.6 * 0.8;
    cases(3).duty_factor    = 1.0;
    cases(3).lam_nm         = 450;
    cases(3).d33_pmV        = 7.0;
    cases(3).n_pump         = 2.05;
    cases(3).n_shg          = 2.05;
    cases(3).overlap_eta    = 0.0152;
    cases(3).width_um       = 0.300;
    cases(3).height_um      = 0.335;
    cases(3).OCE            = 0.15;
    cases(3).uv_loss_mode   = 'multiplier';
    cases(3).uv_loss_val    = 10;
    cases(3).dk_center      = 0;

    % --- CASE 4 ---
    cases(4).name           = 'Input Coupling Improvement (60% -> 80%)';
    cases(4).Pp_avg_mW      = 600 * 0.8 * 0.8;
    cases(4).duty_factor    = 1.0;
    cases(4).lam_nm         = 450;
    cases(4).d33_pmV        = 7.0;
    cases(4).n_pump         = 2.05;
    cases(4).n_shg          = 2.05;
    cases(4).overlap_eta    = 0.0152;
    cases(4).width_um       = 0.300;
    cases(4).height_um      = 0.335;
    cases(4).OCE            = 0.15;
    cases(4).uv_loss_mode   = 'multiplier';
    cases(4).uv_loss_val    = 10;
    cases(4).dk_center      = 0;

    % --- CASE 5 ---
    cases(5).name           = 'Input Facet ARC';
    cases(5).Pp_avg_mW      = 600 * 0.8 * 1.0;
    cases(5).duty_factor    = 1.0;
    cases(5).lam_nm         = 450;
    cases(5).d33_pmV        = 7.0;
    cases(5).n_pump         = 2.05;
    cases(5).n_shg          = 2.05;
    cases(5).overlap_eta    = 0.0152;
    cases(5).width_um       = 0.300;
    cases(5).height_um      = 0.335;
    cases(5).OCE            = 0.15;
    cases(5).uv_loss_mode   = 'multiplier';
    cases(5).uv_loss_val    = 10;
    cases(5).dk_center      = 0;

    % --- CASE 6 ---
    cases(6).name           = 'Pulsing (2x)';
    cases(6).Pp_avg_mW      = 1200 * 0.8 * 1.0;
    cases(6).duty_factor    = 1.0;
    cases(6).lam_nm         = 450;
    cases(6).d33_pmV        = 7.0;
    cases(6).n_pump         = 2.05;
    cases(6).n_shg          = 2.05;
    cases(6).overlap_eta    = 0.0152;
    cases(6).width_um       = 0.300;
    cases(6).height_um      = 0.335;
    cases(6).OCE            = 0.15;
    cases(6).uv_loss_mode   = 'multiplier';
    cases(6).uv_loss_val    = 10;
    cases(6).dk_center      = 0;

    % --- CASE 7 ---
    cases(7).name           = 'Output Coupler';
    cases(7).Pp_avg_mW      = 1200 * 0.8 * 1.0;
    cases(7).duty_factor    = 1.0;
    cases(7).lam_nm         = 450;
    cases(7).d33_pmV        = 7.0;
    cases(7).n_pump         = 2.05;
    cases(7).n_shg          = 2.05;
    cases(7).overlap_eta    = 0.0152;
    cases(7).width_um       = 0.300;
    cases(7).height_um      = 0.335;
    cases(7).OCE            = 0.85;
    cases(7).uv_loss_mode   = 'multiplier';
    cases(7).uv_loss_val    = 10;
    cases(7).dk_center      = 0;

    %% --- 3. LOSS SWEEP CONFIGURATION ---
    blue_losses = [50, 35, 25, 20, 10, 5];   % pump loss values (dB/cm)
    num_losses  = length(blue_losses);
    num_cases   = length(cases);
    num_lengths = length(sweep_vals);

    %% --- 4. EXECUTE SIMULATION ---
    master_guided = zeros(num_lengths, num_losses, num_cases);
    master_scat   = zeros(num_lengths, num_losses, num_cases);

    fprintf('SHG Design Suite v3.0\n');
    fprintf('=====================\n');

    for k = 1:num_cases
        fprintf('\nCase %d/%d: %s\n', k, num_cases, cases(k).name);

        cfg.sim  = sim;
        cfg.phys = cases(k);
        cfg.phys.Pp_peak_mW = cfg.phys.Pp_avg_mW / cfg.phys.duty_factor;

        % Compute g and print overlap info for verification against Wei
        cfg.phys.g = calculate_g(cfg);
        A_eff_m2   = (cfg.phys.width_um * cfg.phys.height_um) * 1e-12;
        ol_factor  = cfg.phys.overlap_eta / sqrt(A_eff_m2);
        fprintf('  overlap_eta = %.4f  (%.2f%%)    ol_factor = %.4ge6 /m    g = %.4e /sqrt(W)/um\n', ...
                cfg.phys.overlap_eta, cfg.phys.overlap_eta * 100, ol_factor / 1e6, cfg.phys.g);

        % Resolve UV loss for each loss level
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

        % Sweep over loss levels and lengths
        for c = 1:num_losses
            case_cfg = cfg;
            case_cfg.phys.a0_dBcm      = loss_cases(c).a0;
            case_cfg.phys.a3_scat_dBcm = loss_cases(c).a3_scat;
            case_cfg.phys.a3_abs_dBcm  = loss_cases(c).a3_abs;

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

    %% --- 6. INTERACTIVE TARGET EVALUATION ---
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

    % Interpolate guided SHG at target length for all cases and loss levels
    matrix_guided = zeros(num_cases, num_losses);
    matrix_scat   = zeros(num_cases, num_losses);
    for k = 1:num_cases
        for c = 1:num_losses
            matrix_guided(k, c) = interp1(sweep_vals, master_guided(:, c, k), target_L, 'linear', 'extrap');
            matrix_scat(k, c)   = interp1(sweep_vals, master_scat(:, c, k),   target_L, 'linear', 'extrap');
        end
    end

    % Print guided SHG table
    fprintf('\n>>> GUIDED SHG POWER  P_guided  (at L = %.0f um) <<<\n', target_L);
    print_table(cases, blue_losses, matrix_guided, num_cases, num_losses);

    % Print scattered SHG table
    fprintf('\n>>> SCATTERED SHG POWER  P_scat  (at L = %.0f um) <<<\n', target_L);
    print_table(cases, blue_losses, matrix_scat, num_cases, num_losses);

    % Print total (guided + scattered) SHG table
    fprintf('\n>>> TOTAL SHG POWER  P_guided + P_scat  (at L = %.0f um) <<<\n', target_L);
    print_table(cases, blue_losses, matrix_guided + matrix_scat, num_cases, num_losses);

    % Optional CSV export
    csv_choice = input('Export guided SHG matrix to CSV? (y/n) [n]: ', 's');
    if strcmpi(csv_choice, 'y') || strcmpi(csv_choice, 'yes')
        filename = sprintf('SHG_v3_0_Evaluation_%.0fum.csv', target_L);
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
%  g = sqrt( 2*w1^2 * deff^2 * ol_factor^2 / (eps0 * c^3 * n1^2 * n2) )
%  ol_factor [1/m] = overlap_eta [dimensionless] / sqrt(A_eff [m^2])
%  Output g is in 1/(sqrt(W)*um) for use with RK4 engine (z in um).
% =========================================================================
function g = calculate_g(cfg)
    eps0   = 8.854e-12;
    c      = 2.998e8;
    lam1   = cfg.phys.lam_nm * 1e-9;
    omega1 = 2 * pi * c / lam1;
    deff   = cfg.phys.d33_pmV * 1e-12;
    n1     = cfg.phys.n_pump;
    n2     = cfg.phys.n_shg;
    A_eff  = (cfg.phys.width_um * cfg.phys.height_um) * 1e-12;   % m^2

    ol_factor = cfg.phys.overlap_eta / sqrt(A_eff);               % 1/m
    yeta      = 2 * omega1^2 / (eps0 * c^3 * n1^2 * n2) * deff^2 * ol_factor^2;
    g_SI      = sqrt(yeta);                                        % 1/(sqrt(W)*m)
    g         = g_SI * 1e-6;                                       % 1/(sqrt(W)*um)
end


%% =========================================================================
%  CORE SOLVER WRAPPER
% =========================================================================
function [Pg, Ps] = run_core_simulation(cfg)
    Pp = cfg.phys.Pp_peak_mW * 1e-3;                              % W (peak)
    a0 = (cfg.phys.a0_dBcm      / 4.3429) * 1e-4;               % Np/um
    as = (cfg.phys.a3_scat_dBcm / 4.3429) * 1e-4;               % Np/um
    aa = (cfg.phys.a3_abs_dBcm  / 4.3429) * 1e-4;               % Np/um

    [Pg, Ps] = rk4_engine(cfg.phys.L_um, cfg.sim.dz, cfg.phys.g, ...
                           a0, as + aa, as, cfg.phys.dk_center, Pp);
end


%% =========================================================================
%  RK4 ENGINE
%  State: Y = [Re(A1), Im(A1), Re(A3), Im(A3), Ps]
%  A1 = pump field amplitude (sqrt(W)), A3 = SH field amplitude (sqrt(W))
%  Ps = accumulated scattered SH power (W)
%  All lengths in um; g in 1/(sqrt(W)*um); losses in Np/um.
% =========================================================================
function [Pg_end, Ps_end] = rk4_engine(L, dz, g, a0, a3t, a3s, dk, Pp)
    z_vec = 0:dz:L;
    Y     = [sqrt(Pp); 0; 0; 0; 0];

    for i = 1:(length(z_vec) - 1)
        z  = z_vec(i);
        k1 = shg_derivs(z,       Y,          a0, a3t, a3s, g, dk);
        k2 = shg_derivs(z+dz/2, Y+k1*dz/2,  a0, a3t, a3s, g, dk);
        k3 = shg_derivs(z+dz/2, Y+k2*dz/2,  a0, a3t, a3s, g, dk);
        k4 = shg_derivs(z+dz,   Y+k3*dz,    a0, a3t, a3s, g, dk);
        Y  = Y + (dz/6) * (k1 + 2*k2 + 2*k3 + k4);
    end

    Pg_end = Y(3)^2 + Y(4)^2;   % guided SH power (W)
    Ps_end = Y(5);               % accumulated scattered SH power (W)
end


function dY = shg_derivs(z, Y, a0, a3t, a3s, g, dk)
    A1 = Y(1) + 1i*Y(2);
    A3 = Y(3) + 1i*Y(4);

    dA1 = -0.5*a0*A1;                                      % undepleted pump (no back-conversion)
    dA3 = -0.5*a3t*A3 + 1i*g*(A1^2)*exp(-1i*dk*z);        % Wei's phase sign convention
    dPs = a3s * abs(A3)^2;                                  % scattered power accumulation

    dY = [real(dA1); imag(dA1); real(dA3); imag(dA3); dPs];
end


%% =========================================================================
%  PER-CASE PLOT
% =========================================================================
function generate_case_plots(case_struct, sweep_vals, guided, scat, blue_losses, param_label)
    num_losses     = length(blue_losses);
    colors         = lines(num_losses);
    h_guided       = gobjects(num_losses, 1);
    legend_labels  = cell(num_losses, 1);

    for c = 1:num_losses
        legend_labels{c} = sprintf('\\alpha_0=%d dB/cm', blue_losses(c));
    end

    fig_title = sprintf('Case: %s', case_struct.name);
    sub_text  = sprintf('P_{avg}=%.1fmW, DF=%.2g (P_{peak}=%.1fmW), UV loss: %s x%g', ...
                case_struct.Pp_avg_mW, case_struct.duty_factor, ...
                case_struct.Pp_avg_mW / case_struct.duty_factor, ...
                case_struct.uv_loss_mode, case_struct.uv_loss_val);

    figure('Color','w','Name',fig_title,'Position',[150, 150, 870, 680]);

    % Subplot 1: Absolute time-averaged power
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

    % Subplot 2: Normalized efficiency
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
%  X-axis: pump loss (dB/cm).  One curve per case.
%  Solid = guided SHG,  Dashed = scattered SHG.
% =========================================================================
function generate_loss_axis_plot(cases, master_guided, master_scat, sweep_vals, blue_losses, target_L_um)
    num_cases  = length(cases);
    num_losses = length(blue_losses);
    colors     = lines(num_cases);

    % Interpolate all cases at target length
    guided_at_L = zeros(num_cases, num_losses);
    scat_at_L   = zeros(num_cases, num_losses);
    for k = 1:num_cases
        for c = 1:num_losses
            guided_at_L(k, c) = interp1(sweep_vals, master_guided(:, c, k), target_L_um, 'linear', 'extrap');
            scat_at_L(k, c)   = interp1(sweep_vals, master_scat(:, c, k),   target_L_um, 'linear', 'extrap');
        end
    end

    % Sort loss axis ascending for clean left-to-right plot
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

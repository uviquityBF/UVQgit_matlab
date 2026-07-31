function SHG_Design_Suite_v2_2()
    close all; clear; clc;

    %% --- 1. GLOBAL SIMULATION CONFIGURATION ---
    sim.use_linewidth  = false;
    sim.do_validation  = false;
    sim.dz             = 0.25;
    
    % --- SWEEP SETUP (Sweep Dimension 1) ---
    study_type = 'length'; 
    sweep_vals = linspace(100, 3000, 50); 
    param_label = 'Waveguide Length [\mum]';

    %% --- 2. DEFINE THE INDEPENDENT INPUT CASES (Requirement 3) ---
    % Add or modify cases here. Each case can have unique power, geometry, or duty cycles.
    
    % --- CASE 1 Configuration ---
    cases(1).name = 'Benchtop RunF (Avo#5)';
    cases(1).Pp_avg_mW          = 10 * 0.05 * 0.8;
    cases(1).duty_factor        = 4e-5;         % duty factor square pulse (Requirement 1)
    cases(1).lam_nm             = 450;
    cases(1).d33_pmV            = 7.0;
    cases(1).n_pump             = 2.15;
    cases(1).n_shg              = 2.6;
    cases(1).overlap_eta        = 0.035;
    cases(1).width_um           = 0.300;
    cases(1).height_um          = 0.335;
    cases(1).OCE                = 0.15;
    cases(1).uv_loss_mode       = 'multiplier'; % 'multiplier' or 'constant' (Requirement 2)
    cases(1).uv_loss_val        = 10;           % 10x multiplier OR flat dB/cm value
    cases(1).dk_center          = 0;

    % --- CASE 2 Configuration ---
    cases(2).name = 'Gen 0 Pkg';
    cases(2).Pp_avg_mW          = 200 * 0.6 * 0.8;
    cases(2).duty_factor        = 1.0;         % duty factor square pulse (Requirement 1)
    cases(2).lam_nm             = 450;
    cases(2).d33_pmV            = 7.0;
    cases(2).n_pump             = 2.15;
    cases(2).n_shg              = 2.6;
    cases(2).overlap_eta        = 0.035;
    cases(2).width_um           = 0.300;
    cases(2).height_um          = 0.335;
    cases(2).OCE                = 0.15;
    cases(2).uv_loss_mode       = 'multiplier'; % 'multiplier' or 'constant' (Requirement 2)
    cases(2).uv_loss_val        = 10;           % 10x multiplier OR flat dB/cm value
    cases(2).dk_center          = 0;

    % --- CASE 3 Configuration ---
    cases(3).name = 'Pump Laser: 600mW';
    cases(3).Pp_avg_mW          = 600 * 0.6 * 0.8;
    cases(3).duty_factor        = 1.0;         % duty factor square pulse (Requirement 1)
    cases(3).lam_nm             = 450;
    cases(3).d33_pmV            = 7.0;
    cases(3).n_pump             = 2.15;
    cases(3).n_shg              = 2.6;
    cases(3).overlap_eta        = 0.035;
    cases(3).width_um           = 0.300;
    cases(3).height_um          = 0.335;
    cases(3).OCE                = 0.15;
    cases(3).uv_loss_mode       = 'multiplier'; % 'multiplier' or 'constant' (Requirement 2)
    cases(3).uv_loss_val        = 10;           % 10x multiplier OR flat dB/cm value
    cases(3).dk_center          = 0;

        % --- CASE 4 Configuration ---
    cases(4).name = 'Input Coupling Improvement 60% -> 80%)';
    cases(4).Pp_avg_mW          = 600 * 0.8 * 0.8;
    cases(4).duty_factor        = 1.0;         % duty factor square pulse (Requirement 1)
    cases(4).lam_nm             = 450;
    cases(4).d33_pmV            = 7.0;
    cases(4).n_pump             = 2.15;
    cases(4).n_shg              = 2.6;
    cases(4).overlap_eta        = 0.035;
    cases(4).width_um           = 0.300;
    cases(4).height_um          = 0.335;
    cases(4).OCE                = 0.15;
    cases(4).uv_loss_mode       = 'multiplier'; % 'multiplier' or 'constant' (Requirement 2)
    cases(4).uv_loss_val        = 10;           % 10x multiplier OR flat dB/cm value
    cases(4).dk_center          = 0;

        % --- CASE 5 Configuration ---
    cases(5).name = 'Input Facet ARC';
    cases(5).Pp_avg_mW          = 600 * 0.8 * 1.0;
    cases(5).duty_factor        = 1.0;         % duty factor square pulse (Requirement 1)
    cases(5).lam_nm             = 450;
    cases(5).d33_pmV            = 7.0;
    cases(5).n_pump             = 2.15;
    cases(5).n_shg              = 2.6;
    cases(5).overlap_eta        = 0.035;
    cases(5).width_um           = 0.300;
    cases(5).height_um          = 0.335;
    cases(5).OCE                = 0.15;
    cases(5).uv_loss_mode       = 'multiplier'; % 'multiplier' or 'constant' (Requirement 2)
    cases(5).uv_loss_val        = 10;           % 10x multiplier OR flat dB/cm value
    cases(5).dk_center          = 0;    
   
        % --- CASE 6 Configuration ---
    cases(6).name = 'Pulsing (2x)';
    cases(6).Pp_avg_mW          = 1200 * 0.8 * 1.0;
    cases(6).duty_factor        = 1.0;         % duty factor square pulse (Requirement 1)
    cases(6).lam_nm             = 450;
    cases(6).d33_pmV            = 7.0;
    cases(6).n_pump             = 2.15;
    cases(6).n_shg              = 2.6;
    cases(6).overlap_eta        = 0.035;
    cases(6).width_um           = 0.300;
    cases(6).height_um          = 0.335;
    cases(6).OCE                = 0.15;
    cases(6).uv_loss_mode       = 'multiplier'; % 'multiplier' or 'constant' (Requirement 2)
    cases(6).uv_loss_val        = 10;           % 10x multiplier OR flat dB/cm value
    cases(6).dk_center          = 0;  
   
         % --- CASE 7 Configuration ---
    cases(7).name = 'Output Coupler';
    cases(7).Pp_avg_mW          = 1200 * 0.8 * 1.0;
    cases(7).duty_factor        = 1.0;         % duty factor square pulse (Requirement 1)
    cases(7).lam_nm             = 450;
    cases(7).d33_pmV            = 7.0;
    cases(7).n_pump             = 2.15;
    cases(7).n_shg              = 2.6;
    cases(7).overlap_eta        = 0.035;
    cases(7).width_um           = 0.300;
    cases(7).height_um          = 0.335;
    cases(7).OCE                = 0.85;
    cases(7).uv_loss_mode       = 'multiplier'; % 'multiplier' or 'constant' (Requirement 2)
    cases(7).uv_loss_val        = 10;           % 10x multiplier OR flat dB/cm value
    cases(7).dk_center          = 0;  
    
    
    % --- MULTI-CASE LOSS SWEEP CONFIGURATION (Sweep Dimension 2) ---
    % Baseline blue pump loss values
    blue_losses = [50, 35, 25, 20, 10, 5]; 
    num_losses = length(blue_losses);
    num_cases = length(cases);
    num_lengths = length(sweep_vals);

    %% --- 3. EXECUTE MASTER SIMULATION ENGINE ---
    % 3D Matrix Allocation to store data for the post-processing prompt
    % Dimensions: [Length_Idx, Loss_Idx, Case_Idx]
    master_guided = zeros(num_lengths, num_losses, num_cases);
    master_scat   = zeros(num_lengths, num_losses, num_cases);

    for k = 1:num_cases
        fprintf('\nRunning Master Case %d/%d: %s...\n', k, num_cases, cases(k).name);
        
        % Build local structural configuration template
        cfg.sim = sim;
        cfg.phys = cases(k);
        
        % PULSED POWER CONVERSION (Requirement 1)
        % Calculate peak instantaneous power from average input power
        cfg.phys.Pp_peak_mW = cfg.phys.Pp_avg_mW / cfg.phys.duty_factor;
        
        % Calculate core coupling constants
        cfg.phys.g_numeric = calculate_g_numeric(cfg);
        cfg.phys.g_analytic = calculate_g_analytic(cfg);
        
        % Generate specific loss profile structures for this case
        loss_cases = struct('a0', {}, 'a3_scat', {}, 'a3_abs', {});
        for c = 1:num_losses
            loss_cases(c).a0 = blue_losses(c);
            
            % Resolve UV Scattering Loss via choice criteria (Requirement 2)
            if strcmp(cfg.phys.uv_loss_mode, 'multiplier')
                loss_cases(c).a3_scat = blue_losses(c) * cfg.phys.uv_loss_val;
            else
                loss_cases(c).a3_scat = cfg.phys.uv_loss_val;
            end
            loss_cases(c).a3_abs = 0; % Internal absorption placeholder
        end
        
        % Run the baseline sweep pipeline
        for c = 1:num_losses
            case_cfg = cfg;
            case_cfg.phys.a0_dBcm      = loss_cases(c).a0;
            case_cfg.phys.a3_scat_dBcm = loss_cases(c).a3_scat;
            case_cfg.phys.a3_abs_dBcm  = loss_cases(c).a3_abs;
            
            for i = 1:num_lengths
                if strcmp(study_type, 'length'), case_cfg.phys.L_um = sweep_vals(i); end
                
                % Call integrator with the elevated instantaneous peak power
                [Pg_peak, Ps_peak] = run_core_simulation(case_cfg);
                
                % BACK-END TIME AVERAGING MATRICES (Requirement 1)
                % Convert peak output power back to time-averaged limits
                master_guided(i, c, k) = Pg_peak * cfg.phys.duty_factor * cfg.phys.OCE;
                master_scat(i, c, k)   = Ps_peak * cfg.phys.duty_factor * cfg.phys.OCE;
            end
        end
        
        %% --- 4. GENERATE INDEPENDENT CASE GRAPH (Requirement 3) ---
        generate_case_plots(cases(k), sweep_vals, master_guided(:,:,k), master_scat(:,:,k), blue_losses, param_label);
    end

    %% --- 5. INTERACTIVE TARGET EVALUATION PROMPT (Requirement 4) ---
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
        fprintf('Invalid or out-of-bounds selection. Defaulting to 2000 um.\n');
        target_L = 2000;
    end
    
    % Extract interpolated cross-sections across all scenarios
    matrix_SHG = zeros(num_cases, num_losses);
    for k = 1:num_cases
        for c = 1:num_losses
            matrix_SHG(k, c) = interp1(sweep_vals, master_guided(:, c, k), target_L, 'linear', 'extrap');
        end
    end
    
    % Print clean formatted tabular results to standard terminal
    fprintf('\n>>> CRITICAL EVALUATION MATRIX: TIME-AVERAGED SHG POWER P_{SHG} (at L = %.1f um) <<<\n', target_L);
    fprintf('%-36s', 'Input Configuration Case Name');
    for c = 1:num_losses
        fprintf(' | a0=%2ddBcm', blue_losses(c));
    end
    fprintf('\n%s\n', repmat('-', 38 + num_losses * 13, 1));
    
    for k = 1:num_cases
        fprintf('%-36s', sprintf('%d: %s', k, cases(k).name));
        for c = 1:num_losses
            fprintf(' | %10.3f uW', matrix_SHG(k, c) * 1e6);
        end
        fprintf('\n');
    end
    fprintf('%s\n', repmat('-', 38 + num_losses * 13, 1));
    
    % Offer automated file tracking dump out to csv
    csv_choice = input('Would you like to export this matrix data to a CSV file? (y/n) [n]: ', 's');
    if strcmpi(csv_choice, 'y') || strcmpi(csv_choice, 'yes')
        filename = sprintf('SHG_Evaluation_Matrix_%.0fum.csv', target_L);
        fid = fopen(filename, 'w');
        
        % Write CSV Header lines
        fprintf(fid, 'Case ID, Case Name');
        for c = 1:num_losses
            fprintf(fid, ',Alpha0_%d_dBcm_uW', blue_losses(c));
        end
        fprintf(fid, '\n');
        
        % Write rows out cleanly
        for k = 1:num_cases
            fprintf(fid, '%d,"%s"', k, cases(k).name);
            for c = 1:num_losses
                fprintf(fid, ',%.6f', matrix_SHG(k, c) * 1e6);
            end
            fprintf(fid, '\n');
        end
        fclose(fid);
        fprintf('SUCCESS: Saved summary matrix cleanly to local spreadsheet file: "%s"\n', filename);
    end
end

%% --- SUB-PLOTTER UTILITY MODULE ---
function generate_case_plots(case_struct, sweep_vals, guided, scat, blue_losses, param_label)
    num_losses = length(blue_losses);
    colors = lines(num_losses);
    
    legend_labels = cell(num_losses * 2, 1);
    for c = 1:num_losses
        legend_labels{2*c-1} = sprintf('\\alpha_0=%d (Guided)', blue_losses(c));
        legend_labels{2*c}   = sprintf('\\alpha_0=%d (Scattered)', blue_losses(c));
    end
    
    % Compile dynamic descriptive titles including Duty Factor tracking markers (Requirement 1)
    fig_title = sprintf('Case Sweep: %s', case_struct.name);
    sub_text = sprintf('Pump P_{avg}=%.1fmW, Duty Factor=%.2f (Peak P_{inst}=%.1fmW), UV Rule: %s (%g)', ...
        case_struct.Pp_avg_mW, case_struct.duty_factor, case_struct.Pp_avg_mW / case_struct.duty_factor, ...
        case_struct.uv_loss_mode, case_struct.uv_loss_val);

    figure('Color','w','Name', fig_title, 'Position', [150, 150, 850, 700]);
    
    % Subplot 1: Absolute Time-Averaged Power
    subplot(2,1,1); hold on;
    for c = 1:num_losses
        semilogy(sweep_vals, guided(:, c)*1e6, 'Color', colors(c,:), 'LineWidth', 1.8);
        semilogy(sweep_vals, scat(:, c)*1e6, '--', 'Color', colors(c,:), 'LineWidth', 1.2);
    end
    set(gca, 'YScale', 'log'); grid on; ylabel('Time-Avg Power [\muW]');
    title({fig_title; sub_text}, 'FontSize', 10);
    legend(legend_labels, 'Location', 'eastoutside');

    % Subplot 2: Intrinsic Efficiencies
    subplot(2,1,2); hold on;
    % Normalized Efficiency calculated relative to average power baseline
    Pp_avg_W = case_struct.Pp_avg_mW * 1e-3;
    for c = 1:num_losses
        Eff_g = (guided(:, c) ./ Pp_avg_W^2) * 100;
        Eff_s = (scat(:, c) ./ Pp_avg_W^2) * 100;
        semilogy(sweep_vals, Eff_g, 'Color', colors(c,:), 'LineWidth', 1.8);
        semilogy(sweep_vals, Eff_s, '--', 'Color', colors(c,:), 'LineWidth', 1.2);
    end
    set(gca, 'YScale', 'log'); grid on; ylabel('Avg Efficiency [%/W]'); xlabel(param_label);
    title('Time-Averaged Intrinsic System Conversion Efficiency', 'FontSize', 10);
    legend(legend_labels, 'Location', 'eastoutside');
end

%% --- NUMERIC INTERFERENCE G-COEFFICIENT ---
function g = calculate_g_numeric(cfg)
    eps0 = 8.854e-12;
    c = 2.998e8;
    lam1 = cfg.phys.lam_nm * 1e-9;
    w1 = 2 * pi * c / lam1;
    d_eff = cfg.phys.d33_pmV * 1e-12;
    A_eff_m2 = (cfg.phys.width_um * cfg.phys.height_um) * 1e-12;
    n1 = cfg.phys.n_pump;
    n2 = cfg.phys.n_shg;
    
    g_SI = (w1 * d_eff / c) * sqrt(2 / (eps0 * n1^2 * n2 * c * A_eff_m2)) * cfg.phys.overlap_eta;
    g = g_SI * 1e-6; 
end

%% --- ANALYTIC WEI G-COEFFICIENT ---
function g = calculate_g_analytic(cfg)
    eps0 = 8.854e-12;
    c = 2.998e8;
    lam2 = (cfg.phys.lam_nm / 2) * 1e-9; 
    w2 = 2 * pi * c / lam2;              
    d_eff = cfg.phys.d33_pmV * 1e-12;    
    A_eff_m2 = (cfg.phys.width_um * cfg.phys.height_um) * 1e-12; 
   
    prefactor = (0.5 * w2^2 * d_eff^2) / (eps0 * cfg.phys.n_pump^2 * cfg.phys.n_shg * c^3);
    g_SI = sqrt(cfg.phys.overlap_eta^2 / A_eff_m2 * prefactor);
    g = g_SI * 1e-6; 
end

%% --- CORE SOLVER WRAPPER ---
function [Pg, Ps] = run_core_simulation(cfg)
    % Core integration engine reads the calculated high-peak instantaneous power
    Pp = cfg.phys.Pp_peak_mW * 1e-3;
    a0 = (cfg.phys.a0_dBcm / 4.3429) * 1e-4;
    as = (cfg.phys.a3_scat_dBcm / 4.3429) * 1e-4;
    aa = (cfg.phys.a3_abs_dBcm / 4.3429) * 1e-4;
    
    [Pg, Ps] = rk4_engine(cfg.phys.L_um, cfg.sim.dz, cfg.phys.g_numeric, a0, as+aa, as, cfg.phys.dk_center, Pp);
end

%% --- THE RK4 ENGINE ---
function [Pg_end, Ps_end] = rk4_engine(L, dz, g, a0, a3t, a3s, dk, Pp)
    z_vec = 0:dz:L;
    Y = [sqrt(Pp); 0; 0; 0; 0]; 
    
    for i = 1:(length(z_vec)-1)
        z = z_vec(i);
        
        k1 = shg_derivs(z,        Y,           a0, a3t, a3s, g, dk);
        k2 = shg_derivs(z + dz/2, Y + k1*dz/2, a0, a3t, a3s, g, dk);
        k3 = shg_derivs(z + dz/2, Y + k2*dz/2, a0, a3t, a3s, g, dk);
        k4 = shg_derivs(z + dz,   Y + k3*dz,   a0, a3t, a3s, g, dk);
        
        Y = Y + (dz/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
    Pg_end = Y(3)^2 + Y(4)^2;
    Ps_end = Y(5);
end

function dY = shg_derivs(z, Y, a0, a3t, a3s, g, dk)
    A1 = Y(1) + 1i*Y(2);
    A3 = Y(3) + 1i*Y(4);
    
    dA1 = -0.5*a0*A1;
    dA3 = -0.5*a3t*A3 + 1i*g*(A1^2)*exp(1i*dk*z);
    dPs = a3s * (abs(A3)^2);
    
    dY = [real(dA1); imag(dA1); real(dA3); imag(dA3); dPs];
end
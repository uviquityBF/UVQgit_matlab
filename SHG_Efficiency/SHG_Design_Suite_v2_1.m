function SHG_Design_Suite_v2_1()
    close all; clear;

    %% --- 1. CONFIGURATION ---
    cfg.sim.use_linewidth  = false;
    cfg.sim.do_validation  = false  % Will validate all cases simultaneously
    cfg.sim.dz             = 0.25;
    
    % --- PHYSICAL INPUTS (Constituent Parts) ---
    cfg.phys.Pp_mW         = 1200*0.8*1.0;
    cfg.phys.lam_nm        = 450;
    cfg.phys.d33_pmV       = 7;      % AlN nonlinear coefficient
    cfg.phys.n_pump        = 2.15;     % n_eff at 450nm
    cfg.phys.n_shg         = 2.6;      % n_eff at 225nm
    cfg.phys.overlap_eta   = 0.035;     % TM00-TM04 Overlap Integral
    cfg.phys.width_um      = 0.300;
    cfg.phys.height_um     = 0.335;
    cfg.phys.disp_slope    = 0.00895;
    cfg.phys.dk_center     = 0;        % Assume perfect phase matching at center
    cfg.phys.OCE           = 0.85;     %output coupling efficiency

    % --- MULTI-CASE LOSS CONFIGURATION (Sweep Dimension 2) ---
    % Define the 5 discrete loss configurations requested
    loss_cases(1).a0 = 5;   loss_cases(1).a3_scat = 50; loss_cases(1).a3_abs = 0;
    loss_cases(2).a0 = 10;  loss_cases(2).a3_scat = 100; loss_cases(2).a3_abs = 0;
    loss_cases(3).a0 = 20;  loss_cases(3).a3_scat = 200; loss_cases(3).a3_abs = 0;
    loss_cases(4).a0 = 25;  loss_cases(4).a3_scat = 250; loss_cases(4).a3_abs = 0;
    loss_cases(5).a0 = 35;  loss_cases(5).a3_scat = 350; loss_cases(5).a3_abs = 0;
    loss_cases(6).a0 = 50;  loss_cases(6).a3_scat = 500; loss_cases(6).a3_abs = 0;
    
    % --- SWEEP SETUP (Sweep Dimension 1) ---
    study_type = 'length'; 
    sweep_vals = linspace(100, 3000, 50); 
    param_label = 'Waveguide Length [\mum]';

    %% --- 2. CALCULATE PHYSICAL G-COEFFICIENTS ---
    cfg.phys.g_numeric = calculate_g_numeric(cfg); 
    cfg.phys.g_analytic = calculate_g_analytic(cfg);
    
    fprintf('Calculated Numeric g:  %.4e [um^-1 * W^-0.5]\n', cfg.phys.g_numeric);
    fprintf('Calculated Analytic g: %.4e [um^-1]\n', cfg.phys.g_analytic);

    %% --- 3. EXECUTE NESTED SWEEPS ---
    num_losses = length(loss_cases);
    num_lengths = length(sweep_vals);
    
    % Allocation matrices: Rows = Length, Columns = Loss Case
    results_guided = zeros(num_lengths, num_losses);
    results_scat   = zeros(num_lengths, num_losses);
    
    for c = 1:num_losses
        % Inject current case losses into baseline configuration
        case_cfg = cfg;
        case_cfg.phys.a0_dBcm      = loss_cases(c).a0;
        case_cfg.phys.a3_scat_dBcm = loss_cases(c).a3_scat;
        case_cfg.phys.a3_abs_dBcm  = loss_cases(c).a3_abs;
        
        for i = 1:num_lengths
            if strcmp(study_type, 'length'), case_cfg.phys.L_um = sweep_vals(i); end
            
            [Pg, Ps] = run_core_simulation(case_cfg);
            results_guided(i, c) = Pg.*cfg.phys.OCE;
            results_scat(i, c)   = Ps.*cfg.phys.OCE;
        end
    end
    
    %% --- 4. COLOR PALETTE DEFINITION ---
    % Generate a high-contrast color array for the lines
    colors = lines(num_losses); 
    legend_labels = cell(num_losses * 2, 1);
    for c = 1:num_losses
        legend_labels{2*c-1} = sprintf('\\alpha_0 = %d dB/cm (Guided)', loss_cases(c).a0);
        legend_labels{2*c}   = sprintf('\\alpha_0 = %d dB/cm (Scattered)', loss_cases(c).a0);
    end

    %% --- 5. MAIN PERFORMANCE PLOTS ---
    Pp_W = cfg.phys.Pp_mW * 1e-3;
    
    figure('Color','w','Name','Loss Case Study Suite v2.1','Position', [100, 100, 800, 700]);
    
    % Plot A: Absolute Power (uW)
    subplot(2,1,1);
    hold on;
    for c = 1:num_losses
        semilogy(sweep_vals, results_guided(:, c)*1e6, 'Color', colors(c,:), 'LineWidth', 1.8);
        semilogy(sweep_vals, results_scat(:, c)*1e6, '--', 'Color', colors(c,:), 'LineWidth', 1.2);
    end
    set(gca, 'YScale', 'log');
    grid on; ylabel('Power [\muW]');
    title('Absolute SHG Power Tracking');
    legend(legend_labels, 'Location', 'eastoutside');

    % Plot B: Normalized Efficiency (%/W)
    subplot(2,1,2);
    hold on;
    for c = 1:num_losses
        Eff_guided = (results_guided(:, c) ./ Pp_W^2) * 100;
        Eff_scat   = (results_scat(:, c) ./ Pp_W^2) * 100;
        
        semilogy(sweep_vals, Eff_guided, 'Color', colors(c,:), 'LineWidth', 1.8);
        semilogy(sweep_vals, Eff_scat, '--', 'Color', colors(c,:), 'LineWidth', 1.2);
    end
    set(gca, 'YScale', 'log');
    grid on; ylabel('Normalized Efficiency [%/W]'); xlabel(param_label);
    title('Intrinsic Conversion Efficiency');
    legend(legend_labels, 'Location', 'eastoutside');

    %% --- 6. OPTIONAL VALIDATION PLOT ---
    if cfg.sim.do_validation && strcmp(study_type, 'length')
        run_validation_plot(cfg, loss_cases, sweep_vals, results_guided);
    end
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
    Pp = cfg.phys.Pp_mW * 1e-3;
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

%% --- MULTI-CASE VALIDATION MODE ---
function run_validation_plot(cfg, loss_cases, L_vals, RK4_guided_matrix)
    num_losses = length(loss_cases);
    Pp = cfg.phys.Pp_mW * 1e-3;
    g = cfg.phys.g_analytic; 
    colors = lines(num_losses);

    figure('Name','v2.1 Multi-Case Validation Engine','Color','w','Position', [150, 150, 750, 500]);
    hold on;
    
    % Loop through each case to verify that the RK4 arrays match the analytical formula
    for c = 1:num_losses
        a0 = (loss_cases(c).a0 / 4.3429) * 1e-4;
        a3 = ((loss_cases(c).a3_scat + loss_cases(c).a3_abs) / 4.3429) * 1e-4;
        da = (a3/2 - a0);
        
        DeRate = exp(-a3 .* L_vals) .* ( (exp(da .* L_vals) - 1).^2 ) ./ ( (da .* L_vals).^2 + 1e-20 );
        P_analytic = (g .* L_vals).^2 * (Pp^2) .* DeRate .* cfg.phys.OCE;

        % Plot numerical values as discrete points
        plot(L_vals, RK4_guided_matrix(:, c)*1e6, 'o', 'Color', colors(c,:), 'MarkerSize', 5);
        % Plot analytical values as clean solid continuous lines
        plot(L_vals, P_analytic*1e6, '-', 'Color', colors(c,:), 'LineWidth', 1.5);
    end
    
    grid on; set(gca, 'YScale', 'log');
    xlabel('Length [\mum]'); ylabel('Guided Power [\muW]');
    title('Multi-Case Validation (Dots = RK4, Lines = Wei Analytical)');
    
    % Truncated visual legend tracking the solid lines
    short_labels = cell(num_losses, 1);
    for c = 1:num_losses
        short_labels{c} = sprintf('\\alpha_0 = %d dB/cm', loss_cases(c).a0);
    end
    legend(short_labels, 'Location', 'eastoutside');
end
function SHG_Design_Suite_v2()
    close all; clear;

    %% --- 1. CONFIGURATION ---
    cfg.sim.use_linewidth  = false;
    cfg.sim.do_validation  = true;
    cfg.sim.dz             = 0.25;
    
    % --- PHYSICAL INPUTS (Constituent Parts) -
    cfg.phys.Pp_mW         = 25;
    cfg.phys.lam_nm        = 450;
    cfg.phys.d33_pmV       = 4.7;      % AlN nonlinear coefficient
    cfg.phys.n_pump        = 2.15;     % n_eff at 450nm
    cfg.phys.n_shg         = 2.6;     % n_eff at 225nm
    cfg.phys.overlap_eta   = 0.03;     % TM00-TM04 Overlap Integral
    cfg.phys.width_um      = 0.300;
    cfg.phys.height_um     = 0.335;
    cfg.phys.disp_slope    = 0.00895;
    cfg.phys.dk_center     = 0;        % Assume perfect phase matching at center

    % --- LOSSES ---
    cfg.phys.a0_dBcm       = 30;
    cfg.phys.a3_scat_dBcm  = 300;
    cfg.phys.a3_abs_dBcm   = 0;
    
    % --- SWEEP SETUP ---
    study_type = 'length'; 
    sweep_vals = linspace(100, 3000, 50); 
    param_label = 'Waveguide Length [\mum]';

    %% --- 2. CALCULATE PHYSICAL G-COEFFICIENTS ---
    cfg.phys.g_numeric = calculate_g_numeric(cfg); 
    cfg.phys.g_analytic = calculate_g_analytic(cfg);
    
    fprintf('Calculated Numeric g:  %.4e [um^-1 * W^-0.5]\n', cfg.phys.g_numeric);
    fprintf('Calculated Analytic g: %.4e [um^-1]\n', cfg.phys.g_analytic);

    %% --- 3. EXECUTE SWEEP ---
    results = zeros(length(sweep_vals), 2); % [Guided, Scattered]
    
    for i = 1:length(sweep_vals)
        current_cfg = cfg;
        if strcmp(study_type, 'length'), current_cfg.phys.L_um = sweep_vals(i); end
        
        [Pg, Ps] = run_core_simulation(current_cfg);
        results(i, :) = [Pg, Ps];
    end
    
    %% --- 4. OPTIONAL VALIDATION PLOT ---
    if cfg.sim.do_validation && strcmp(study_type, 'length')
        run_validation_plot(cfg, sweep_vals, results(:,1));
    end
    
    %% --- 5. CONVERT TO NORMALIZED EFFICIENCY (%/W) ---
    Pp_W = cfg.phys.Pp_mW * 1e-3;
    Eff_guided = (results(:,1) ./ Pp_W^2) * 100;
    Eff_scat   = (results(:,2) ./ Pp_W^2) * 100;

    %% --- 6. PLOTTING ---
    figure('Color','w','Name','Normalized Efficiency Study');
    
    % Plot A: Absolute Power (uW)
    subplot(2,1,1);
    semilogy(sweep_vals, results(:,1)*1e6, 'b', 'LineWidth', 1.5); hold on;
    semilogy(sweep_vals, results(:,2)*1e6, 'r--', 'LineWidth', 1.5);
    grid on; ylabel('Power [\muW]'); legend('Guided','Scattered');
    title('Absolute SHG Power');

    % Plot B: Normalized Efficiency (%/W)
    subplot(2,1,2);
    plot(sweep_vals, Eff_guided, 'b', 'LineWidth', 1.5); hold on;
    plot(sweep_vals, Eff_scat, 'r--', 'LineWidth', 1.5);
    grid on; ylabel('Normalized Efficiency [%/W]'); xlabel(param_label);
    title('Intrinsic Conversion Efficiency');
end

%% --- NUMERIC INTERFERENCE G-COEFFICIENT (For RK4 Tracking) ---
function g = calculate_g_numeric(cfg)
    eps0 = 8.854e-12;
    c = 2.998e8;
    
    lam1 = cfg.phys.lam_nm * 1e-9;
    w1 = 2 * pi * c / lam1;
    
    d_eff = cfg.phys.d33_pmV * 1e-12;
    A_eff_m2 = (cfg.phys.width_um * cfg.phys.height_um) * 1e-12;
    
    n1 = cfg.phys.n_pump;
    n2 = cfg.phys.n_shg;
    
    % Core waveguide coupling factor for direct root-power arrays
    g_SI = (w1 * d_eff / c) * sqrt(2 / (eps0 * n1^2 * n2 * c * A_eff_m2)) * cfg.phys.overlap_eta;
    g = g_SI * 1e-6; % m^-1 to um^-1
    
    %g = g*sqrt(8);
end

%% --- ANALYTIC WEI G-COEFFICIENT (For Validation Equation Only) ---
function g = calculate_g_analytic(cfg)
    eps0 = 8.854e-12;
    c = 2.998e8;
    
    lam2 = (cfg.phys.lam_nm / 2) * 1e-9; 
    w2 = 2 * pi * c / lam2;              
    d_eff = cfg.phys.d33_pmV * 1e-12;    
    A_eff_m2 = (cfg.phys.width_um * cfg.phys.height_um) * 1e-12; 
   
    % Traditional plane-wave prefactor carrying implicit power normalization matching Wei requirements
    %     prefactor = (2 * w2^2 * d_eff^2) / (eps0 * cfg.phys.n_pump^2 * cfg.phys.n_shg * c^3);
    % REMOVE THE PEAK-FACTOR OVER-SCALING: 
    % Changing the prefactor from 2 to 0.5 strips the legacy continuous plane-wave 
    % peak envelope fraction, shifting the analytical model into the exact same 
    % RMS root-power space (|A|^2 = Watts) used by your RK4 engine.
    prefactor = (0.5 * w2^2 * d_eff^2) / (eps0 * cfg.phys.n_pump^2 * cfg.phys.n_shg * c^3);
    
    g_SI = sqrt(cfg.phys.overlap_eta^2 / A_eff_m2 * prefactor);
    g = g_SI * 1e-6; % m^-1 to um^-1   
   % g = g*sqrt(8);

end

%% --- CORE SOLVER WRAPPER ---
function [Pg, Ps] = run_core_simulation(cfg)
    Pp = cfg.phys.Pp_mW * 1e-3;
    a0 = (cfg.phys.a0_dBcm / 4.3429) * 1e-4;
    as = (cfg.phys.a3_scat_dBcm / 4.3429) * 1e-4;
    aa = (cfg.phys.a3_abs_dBcm / 4.3429) * 1e-4;
    
    % Always pass the specialized numeric tracking g to the integrator grid
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

%% --- VALIDATION MODE ---
function run_validation_plot(cfg, L_vals, RK4_guided)
    a0 = (cfg.phys.a0_dBcm / 4.3429) * 1e-4;
    a3 = ((cfg.phys.a3_scat_dBcm + cfg.phys.a3_abs_dBcm) / 4.3429) * 1e-4;
    Pp = cfg.phys.Pp_mW * 1e-3;
    
    % Explicitly draw the analytic efficiency coefficient for the continuous equations
    g = cfg.phys.g_analytic; 
    da = (a3/2 - a0);
    
    DeRate = exp(-a3 .* L_vals) .* ( (exp(da .* L_vals) - 1).^2 ) ./ ( (da .* L_vals).^2 + 1e-20 );
    
    % Pure analytical Wei solution utilizing single power scalar dependency
    P_analytic = (g .* L_vals).^2 * (Pp^2) .* DeRate;

    figure('Name','RK4 vs Wei Validation','Color','w');
    semilogy(L_vals, RK4_guided*1e6, 'bo', 'LineWidth', 1.5, 'DisplayName', 'RK4 Numerical'); hold on;
    semilogy(L_vals, P_analytic*1e6, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Wei Analytical');
    grid on; xlabel('Length [\mum]'); ylabel('Guided Power [\muW]');
    legend('Location','best'); title('Solver Validation (Analytical vs Numerical)');
end




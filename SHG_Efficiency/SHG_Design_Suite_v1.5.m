function SHG_Design_Suite_v2()
    close all; clear;

    %% --- 1. CONFIGURATION ---
    cfg.sim.use_linewidth  = false;
    cfg.sim.do_validation  = true;
    cfg.sim.dz             = 0.5;
    
    % --- PHYSICAL INPUTS (Constituent Parts) ---
    cfg.phys.Pp_mW         = 100;
    cfg.phys.lam_nm        = 450;
    cfg.phys.d33_pmV       = 4.7;      % AlN nonlinear coefficient
    cfg.phys.n_pump        = 2.15;     % n_eff at 450nm
    cfg.phys.n_shg         = 2.25;     % n_eff at 225nm
    cfg.phys.overlap_eta   = 0.03;     % TM00-TM04 Overlap Integral
    cfg.phys.width_um      = 0.3365;
    cfg.phys.height_um     = 0.381;
    cfg.phys.disp_slope    = 0.0089;
    cfg.phys.dk_center     = 0;        % Assume perfect phase matching at center

    % --- LOSSES ---
    cfg.phys.a0_dBcm       = 5;
    cfg.phys.a3_scat_dBcm  = 200;
    cfg.phys.a3_abs_dBcm   = 20;
    
    % --- SWEEP SETUP ---
    study_type = 'length'; 
    sweep_vals = linspace(100, 5000, 40); 
    param_label = 'Waveguide Length [\mum]';

    %% --- 2. CALCULATE PHYSICAL G-COEFFICIENT ---
    cfg.phys.g_coeff = calculate_g(cfg); % No longer a hardcoded number!
    fprintf('Calculated g_coeff: %.4e [um^-1]\n', cfg.phys.g_coeff);

    %% --- 3. EXECUTE SWEEP ---
    results = zeros(length(sweep_vals), 2); % [Guided, Scattered]
    
    for i = 1:length(sweep_vals)
        current_cfg = cfg;
        if strcmp(study_type, 'length'), current_cfg.phys.L_um = sweep_vals(i); end
        
        [Pg, Ps] = run_core_simulation(current_cfg);
        results(i, :) = [Pg, Ps];
    end

    %% --- 4. CONVERT TO NORMALIZED EFFICIENCY (%/W) ---
    Pp_W = cfg.phys.Pp_mW * 1e-3;
    % Efficiency = (P_shg / P_pump^2) * 100
    Eff_guided = (results(:,1) ./ Pp_W^2) * 100;
    Eff_scat   = (results(:,2) ./ Pp_W^2) * 100;

    %% --- 5. PLOTTING ---
    figure('Color','w','Name','Normalized Efficiency Study');
    
    % Plot A: Absolute Power (uW)
    subplot(2,1,1);
    plot(sweep_vals, results(:,1)*1e6, 'b', 'LineWidth', 1.5); hold on;
    plot(sweep_vals, results(:,2)*1e6, 'r--', 'LineWidth', 1.5);
    grid on; ylabel('Power [\muW]'); legend('Guided','Scattered');
    title('Absolute SHG Power');

    % Plot B: Normalized Efficiency (%/W)
    subplot(2,1,2);
    plot(sweep_vals, Eff_guided, 'b', 'LineWidth', 1.5); hold on;
    plot(sweep_vals, Eff_scat, 'r--', 'LineWidth', 1.5);
    grid on; ylabel('Normalized Efficiency [%/W]'); xlabel(param_label);
    title('Intrinsic Conversion Efficiency');
end

%% --- CONSTITUENT G-COEFFICIENT CALCULATION ---
function g = calculate_g(cfg)
    % Based on Boyd's Nonlinear Optics (Coupled Wave Equations)
    % g = sqrt( (2 * w2^2 * d_eff^2) / (eps0 * n1^2 * n2 * c^3 * A_eff) ) * overlap
    
    eps0 = 8.854e-12;
    c = 2.998e8;
    lam2 = (cfg.phys.lam_nm / 2) * 1e-9; % SHG wavelength in meters
    w2 = 2 * pi * c / lam2;              % SHG angular frequency
    d_eff = cfg.phys.d33_pmV * 1e-12;    % Convert pm/V to m/V
    
    A_eff_m2 = (cfg.phys.width_um * cfg.phys.height_um) * 1e-12; % um^2 to m^2
    
    % The prefactor term (SI Units)
    prefactor = (2 * w2^2 * d_eff^2) / (eps0 * cfg.phys.n_pump^2 * cfg.phys.n_shg * c^3);
    
    % Combine with overlap and area, then convert back to um^-1 for the integrator
    % g_SI has units of [m^-1 / sqrt(W)]
    g_SI = sqrt(cfg.phys.overlap_eta^2 / A_eff_m2 * prefactor);
    
    g = g_SI * 1e-6; % Convert m^-1 to um^-1
end

%% --- CORE SOLVER WRAPPER ---
function [Pg, Ps] = run_core_simulation(cfg)
    Pp = cfg.phys.Pp_mW * 1e-3;
    a0 = (cfg.phys.a0_dBcm / 4.3429) * 1e-4;
    as = (cfg.phys.a3_scat_dBcm / 4.3429) * 1e-4;
    aa = (cfg.phys.a3_abs_dBcm / 4.3429) * 1e-4;
    
    if cfg.sim.use_linewidth
        % Spectral Overlap logic (Omitted here for brevity, reuse from previous)
        [Pg, Ps] = handle_spectral_overlap(cfg, a0, as, aa, Pp);
    else
        [Pg, Ps] = rk4_engine(cfg.phys.L_um, cfg.sim.dz, cfg.phys.g_coeff, a0, as+aa, as, cfg.phys.dk_center, Pp);
    end
end

%% --- THE RK4 ENGINE ---
function [Pg_end, Ps_end] = rk4_engine(L, dz, g, a0, a3t, a3s, dk, Pp)
    z_vec = 0:dz:L;
    Y = [sqrt(Pp); 0; 0; 0; 0; 0]; % [A1r, A1i, A3r, A3i, Pscat, Phase]
    for i = 1:(length(z_vec)-1)
        k1 = shg_derivs(Y, a0, a3t, a3s, g, dk);
        k2 = shg_derivs(Y + k1*dz/2, a0, a3t, a3s, g, dk);
        k3 = shg_derivs(Y + k2*dz/2, a0, a3t, a3s, g, dk);
        k4 = shg_derivs(Y + k3*dz, a0, a3t, a3s, g, dk);
        Y = Y + (dz/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
    Pg_end = Y(3)^2 + Y(4)^2;
    Ps_end = Y(5);
end

function dY = shg_derivs(Y, a0, a3t, a3s, g, dk)
    A1 = Y(1) + 1i*Y(2);
    A3 = Y(3) + 1i*Y(4);
    dA1 = -0.5*a0*A1;
    dA3 = -0.5*a3t*A3 + 1i*g*(A1^2)*exp(-1i*Y(6));
    dY = [real(dA1); imag(dA1); real(dA3); imag(dA3); a3s*abs(A3)^2; dk];
end

%% --- VALIDATION MODE ---
function run_validation_plot(cfg, L_vals, RK4_guided)
    % Calculate Wei Analytical Solution for comparison
    a0 = (cfg.phys.a0_dBcm / 4.3429) * 1e-4;
    a3 = ((cfg.phys.a3_scat_dBcm + cfg.phys.a3_abs_dBcm) / 4.3429) * 1e-4;
    Pp = cfg.phys.Pp_mW * 1e-3;
    g = cfg.phys.g_coeff;
    
    L_cm = L_vals * 1e-4;
    da = (a3/2 - a0);
    % Wei Formula
    DeRate = exp(-a3*L_vals*1e-4) .* ( (exp(da*L_vals*1e-4) - 1).^2 ) ./ ( (da*L_vals*1e-4).^2 + 1e-20 );
    P_analytic = (g * L_vals).^2 * Pp .* DeRate;

    figure('Name','RK4 vs Wei Validation');
    plot(L_vals, RK4_guided*1e6, 'bo', 'DisplayName', 'RK4 Numerical'); hold on;
    plot(L_vals, P_analytic*1e6, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Wei Analytical');
    grid on; xlabel('Length [\mum]'); ylabel('Guided Power [\muW]');
    legend; title('Solver Validation (Analytical vs Numerical)');
end
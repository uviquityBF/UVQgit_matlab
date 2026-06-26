function SHG_Design_Suite()
    % SHG_Design_Suite: Advanced RK4 Nonlinear Simulator
    % Focus: Absolute Power, Validation, and Multi-parameter Sweeps.
    
    close all; clear;

    %% --- 1. GLOBAL CONFIGURATION ---
    cfg.sim.use_linewidth  = false;      % Speed toggle for spectral overlap
    cfg.sim.do_validation  = true;       % Overlay Analytical vs RK4 for 1st case
    cfg.sim.dz             = 0.5;        % [um] Integrator step
    
    % Default Physical Constants (AlN TM04 example)
    cfg.phys.Pp_mW         = 100;
    cfg.phys.lam_nm        = 450;
    cfg.phys.linewidth_nm  = 0.15;
    cfg.phys.g_coeff       = 1.15e-4;    % Calculated nonlinear coupling [um^-1]
    cfg.phys.disp_slope    = 0.0089;     % Delta-n slope [1/nm]
    cfg.phys.dk_center     = 0;          % Initial detuning [um^-1]

    % Default Loss/Geo (Will be overridden during sweeps)
    cfg.phys.a0_dBcm       = 5;
    cfg.phys.a3_scat_dBcm  = 200;
    cfg.phys.a3_abs_dBcm   = 20;
    cfg.phys.L_um          = 2000;

    %% --- 2. SELECT STUDY TYPE ---
    % Options: 'length', 'pump_loss', 'shg_scat', 'shg_abs'
    study_type = 'length'; 
    
    switch study_type
        case 'length'
            sweep_vals = linspace(100, 5000, 30); % [um]
            param_label = 'Waveguide Length [\mum]';
        case 'pump_loss'
            sweep_vals = linspace(0, 50, 20);      % [dB/cm]
            param_label = 'Pump Loss [dB/cm]';
        case 'shg_scat'
            sweep_vals = logspace(1, 3.5, 25);    % [dB/cm]
            param_label = 'SHG Scatter Rate [dB/cm]';
    end

    %% --- 3. EXECUTE SWEEP ---
    results = zeros(length(sweep_vals), 3); % [Guided, Scattered, Total]
    
    fprintf('Running %s sweep...\n', study_type);
    for i = 1:length(sweep_vals)
        % Update the specific parameter being swept
        current_cfg = cfg;
        if strcmp(study_type, 'length'),    current_cfg.phys.L_um = sweep_vals(i); end
        if strcmp(study_type, 'pump_loss'), current_cfg.phys.a0_dBcm = sweep_vals(i); end
        if strcmp(study_type, 'shg_scat'),  current_cfg.phys.a3_scat_dBcm = sweep_vals(i); end
        
        % Run the Solver
        [Pg, Ps] = run_core_simulation(current_cfg);
        results(i, :) = [Pg, Ps, Pg+Ps];
    end

    %% --- 4. OPTIONAL VALIDATION PLOT ---
    if cfg.sim.do_validation && strcmp(study_type, 'length')
        run_validation_plot(cfg, sweep_vals, results(:,1));
    end

    %% --- 5. KEY FIGURES OF MERIT PLOTS ---
    figure('Color','w','Name','Parametric Study Results');
    subplot(2,1,1);
    plot(sweep_vals, results(:,1)*1e6, 'b-o', 'LineWidth', 1.5); hold on;
    plot(sweep_vals, results(:,2)*1e6, 'r-s', 'LineWidth', 1.5);
    grid on; ylabel('Power [\muW]'); legend('Guided','Scattered');
    title(['Absolute SHG Power vs ', param_label]);

    subplot(2,1,2);
    plot(sweep_vals, results(:,2)./results(:,1), 'k', 'LineWidth', 2);
    grid on; ylabel('Ratio (Scat/Guided)'); xlabel(param_label);
    title('Extraction Efficiency Ratio');
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
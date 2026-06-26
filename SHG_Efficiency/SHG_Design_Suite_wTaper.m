function SHG_Design_Suite_wTaper()
    close all; clear;

    %% --- 1. CONFIGURATION ---
    cfg.sim.use_linewidth  = false;
    cfg.sim.do_validation  = false;
    cfg.sim.dz             = 0.5;
    
    % --- PHYSICAL INPUTS ---
    cfg.phys.Pp_mW         = 100;
    cfg.phys.lam_nm        = 450;
    cfg.phys.d33_pmV       = 4.7;      
    cfg.phys.n_pump        = 2.15;     
    cfg.phys.n_shg         = 2.25;     
    cfg.phys.overlap_eta   = 0.03;     
    cfg.phys.width_um      = 0.3365;  % Nominal width
    cfg.phys.height_um     = 0.381;
    cfg.phys.dk_center     = 0;        

    % --- GEOMETRY / TAPER SETUP ---
    cfg.geo.L_um           = 2000;
    cfg.geo.w_nominal      = 0.3365;
    cfg.geo.w_start        = 0.3365;  % Change these to create a taper
    cfg.geo.w_end          = 0.3365;  % (Keep same for validation)
    cfg.geo.taper_slope    = 0.5;     % dk shift per micron of width change

    % --- LOSSES ---
    cfg.phys.a0_dBcm       = 20;
    cfg.phys.a3_scat_dBcm  = cfg.phys.a0_dBcm * 10;
    cfg.phys.a3_abs_dBcm   = 0;
    
    % --- SWEEP SETUP ---
    study_type = 'length'; 
    sweep_vals = linspace(100, 2000, 40); 
    param_label = 'Waveguide Length [\mum]';

    %% --- 2. PRE-CALCS ---
    cfg.phys.g_coeff = calculate_g(cfg); 
    fprintf('Calculated g_coeff: %.4e [um^-1]\n', cfg.phys.g_coeff);

    %% --- 3. EXECUTE SWEEP ---
    results = zeros(length(sweep_vals), 2); 
    
    for i = 1:length(sweep_vals)
        current_cfg = cfg;
        if strcmp(study_type, 'length'), current_cfg.geo.L_um = sweep_vals(i); end
        
        [Pg, Ps] = run_core_simulation(current_cfg);
        results(i, :) = [Pg, Ps];
    end

    %% --- 4. CONVERT TO NORMALIZED EFFICIENCY ---
    Pp_W = cfg.phys.Pp_mW * 1e-3;
    Eff_guided = (results(:,1) ./ Pp_W^2) * 100;
    Eff_scat   = (results(:,2) ./ Pp_W^2) * 100;

%% --- 5. PLOTTING (Semilogy Version) ---
    fig_title = sprintf('SHG Study: Pp=%.1fmW, g=%.2e, a0=%.1fdB, as=%.1fdB, aa=%.1fdB', ...
                cfg.phys.Pp_mW, cfg.phys.g_coeff, cfg.phys.a0_dBcm, ...
                cfg.phys.a3_scat_dBcm, cfg.phys.a3_abs_dBcm);
            
    figure('Color', 'w', 'Name', 'SHG Results', 'Position', [100 100 800 700]);
    
    L_vec = sweep_vals;
    P_guide = results(:,1) * 1e6; % uW
    P_scat  = results(:,2) * 1e6; % uW
    
    % --- Plot A: Absolute Power (uW) ---
    subplot(2,1,1);
    h1 = semilogy(L_vec, P_guide, 'b', 'LineWidth', 2); hold on;
    h2 = semilogy(L_vec, P_scat, 'r--', 'LineWidth', 2);
    grid on; ylabel('Power [\muW] (log)');
    title({'Absolute SHG Power', fig_title}, 'FontSize', 10);
    legend('Guided SHG', 'Scattered SHG', 'Location', 'southwest');

    % Peak and End Labels
    [maxG, idxG] = max(P_guide);
    [maxS, idxS] = max(P_scat);
    
    % Use a small offset for log-scale text positioning
    text(L_vec(idxG), maxG*1.2, sprintf(' Peak: %.2f\\muW', maxG), 'Color', 'b');
    text(L_vec(end), P_guide(end), sprintf(' End: %.2f\\muW', P_guide(end)), 'Color', 'b', 'HorizontalAlignment', 'left');
    plot(L_vec(idxG), maxG, 'bo', 'MarkerFaceColor', 'b');

    text(L_vec(idxS), maxS*0.8, sprintf(' Peak: %.2f\\muW', maxS), 'Color', 'r', 'VerticalAlignment', 'top');
    text(L_vec(end), P_scat(end), sprintf(' End: %.2f\\muW', P_scat(end)), 'Color', 'r', 'HorizontalAlignment', 'left');
    plot(L_vec(idxS), maxS, 'rs', 'MarkerFaceColor', 'r');
    
    % --- Plot B: Normalized Efficiency (%/W) ---
    subplot(2,1,2);
    semilogy(L_vec, Eff_guided, 'b', 'LineWidth', 2); hold on;
    semilogy(L_vec, Eff_scat, 'r--', 'LineWidth', 2);
    grid on; ylabel('Efficiency [%/W] (log)'); xlabel(param_label);
    title('Normalized Conversion Efficiency');
    
    [maxEG, idxEG] = max(Eff_guided);
    text(L_vec(idxEG), maxEG*1.2, sprintf(' %.3f%%/W', maxEG), 'Color', 'b');
    text(L_vec(end), Eff_guided(end), sprintf(' %.3f', Eff_guided(end)), 'Color', 'b');
end

%% --- PHYSICAL CALCULATIONS ---
function g = calculate_g(cfg)
    eps0 = 8.854e-12;
    c = 2.998e8;
    lam2 = (cfg.phys.lam_nm / 2) * 1e-9; 
    w2 = 2 * pi * c / lam2;              
    d_eff = cfg.phys.d33_pmV * 1e-12;    
    A_eff_m2 = (cfg.phys.width_um * cfg.phys.height_um) * 1e-12; 
    prefactor = (2 * w2^2 * d_eff^2) / (eps0 * cfg.phys.n_pump^2 * cfg.phys.n_shg * c^3);
    g_SI = 2*sqrt(cfg.phys.overlap_eta^2 / A_eff_m2 * prefactor);
    g = g_SI * 1e-6; 
end

function [Pg, Ps] = run_core_simulation(cfg)
    Pp = cfg.phys.Pp_mW * 1e-3;
    a0 = (cfg.phys.a0_dBcm / 4.3429) * 1e-4;
    as = (cfg.phys.a3_scat_dBcm / 4.3429) * 1e-4;
    aa = (cfg.phys.a3_abs_dBcm / 4.3429) * 1e-4;
    
    % Use the engine and pass the whole cfg struct
    [Pg, Ps] = rk4_engine_tapered(cfg.geo.L_um, cfg.sim.dz, cfg, a0, as+aa, as, Pp);
end

function [Pg_end, Ps_end] = rk4_engine_tapered(L, dz, cfg, a0, a3t, a3s, Pp)
    z_vec = 0:dz:L;
    Y = [sqrt(Pp); 0; 0; 0; 0; 0];  
    for i = 1:(length(z_vec)-1)
        z_now = z_vec(i);
        k1 = shg_derivs_linear_taper(z_now,        Y,           a0, a3t, a3s, cfg);
        k2 = shg_derivs_linear_taper(z_now + dz/2, Y + k1*dz/2, a0, a3t, a3s, cfg);
        k3 = shg_derivs_linear_taper(z_now + dz/2, Y + k2*dz/2, a0, a3t, a3s, cfg);
        k4 = shg_derivs_linear_taper(z_now + dz,   Y + k3*dz,   a0, a3t, a3s, cfg);
        Y = Y + (dz/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
    Pg_end = Y(3)^2 + Y(4)^2;
    Ps_end = Y(5);
end

function dY = shg_derivs_linear_taper(z, Y, a0, a3t, a3s, cfg)
    % Linear Taper Width
    L_eff = max(cfg.geo.L_um, 1e-9); % Prevent div by zero
    w_z = cfg.geo.w_start + (cfg.geo.w_end - cfg.geo.w_start) * (z / L_eff);
    
    % Linear Phase Mismatch Approximation
    dw = w_z - cfg.geo.w_nominal;
    dk_local = cfg.phys.dk_center + (dw * cfg.geo.taper_slope);
    
    A1 = Y(1) + 1i*Y(2);
    A3 = Y(3) + 1i*Y(4);
    phi_acc = Y(6); 
    
    dA1 = -0.5 * a0 * A1;
    dA3 = -0.5 * a3t * A3 + 1i * cfg.phys.g_coeff * (A1^2) * exp(-1i * phi_acc);
    
    dY = [real(dA1); imag(dA1); real(dA3); imag(dA3); a3s*abs(A3)^2; dk_local];
end

% function run_validation_plot(cfg, L_vals, RK4_guided)
%     Pp = (cfg.phys.Pp_mW * 1e-3);
%     g = cfg.phys.g_coeff;
%     
%     % IDEAL CASE (No loss, no dk)
%     P_ideal = (g * L_vals).^2 * Pp^2; 
% 
%     figure('Name','The Vacuum Test');
%     plot(L_vals, RK4_guided*1e6, 'bo', 'DisplayName', 'RK4 Numerical'); hold on;
%     plot(L_vals, P_ideal*1e6, 'r-', 'DisplayName', 'Ideal Quadratic');
%     legend; grid on;
% end

function run_validation_plot(cfg, L_vals, RK4_guided)
    % 1. Convert losses to per-micron (matching RK4)
    a0 = (cfg.phys.a0_dBcm / 4.3429) * 1e-4; 
    as = (cfg.phys.a3_scat_dBcm / 4.3429) * 1e-4;
    aa = (cfg.phys.a3_abs_dBcm / 4.3429) * 1e-4;
    a3 = as + aa;
    
    Pp = (cfg.phys.Pp_mW * 1e-3);
    g = cfg.phys.g_coeff;
    
    % 2. Numerical Stability Check
    % Instead of (exp(da*L)-1)^2 / da^2, we use the hyperbolic sinc
    % which is much more stable as da approaches 0.
    da = (a3/2 - a0);
    
    % We use the fact that (exp(x)-1)/x = exp(x/2) * sinh(x/2)/(x/2)
    % This is exactly equivalent to your Wei formula but won't crash at zero.
    term1 = exp(-a3 * L_vals);
    term2 = exp(da * L_vals);
    
% Define the argument x = da * L
    x = da * L_vals;
    
    % Initialize growth_factor
    growth_f = zeros(size(x));
    
    % Use Taylor expansion for stability if x is very small (near 0)
    % otherwise use the standard formula.
    small_idx = abs(x) < 1e-4;
    
    % Standard formula
    growth_f(~small_idx) = ((exp(x(~small_idx)) - 1) ./ x(~small_idx)).^2;
    
    % Taylor expansion: (1 + x/2 + x^2/6 + ...)^2
    growth_f(small_idx) = (1 + x(small_idx)/2 + (x(small_idx).^2)/6).^2;
    
    P_analytic = (g * L_vals).^2 * Pp^2 .* exp(-a3 * L_vals) .* growth_f;
    
    % Re-calculating P_analytic
    % P = (g*L)^2 * Pp^2 * exp(-a3*L) * growth_factor
    P_analytic = (g * L_vals).^2 * Pp^2 .* term1 .* growth_f;

    figure('Name','RK4 vs Wei Validation');
    plot(L_vals, RK4_guided*1e6, 'bo', 'DisplayName', 'RK4 Numerical'); hold on;
    plot(L_vals, P_analytic*1e6, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Wei Analytical');
    grid on; xlabel('Length [\mum]'); ylabel('Guided Power [\muW]');
    legend; title('Validation (Numerical vs. Analytical)');
end
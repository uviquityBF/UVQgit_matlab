function SHG_Wei_Reconciliation()
% Standalone reconciliation script: Brent's SHG model vs Wei's benchmark.
%
% Parameters are set to match Wei's benchmark (SHG_AlN_UVC_CME_benchmark_final.m).
%
% Four methods are computed to isolate each source of disagreement:
%   Method 1 (Wei Analytic):        Wei's closed-form PM_factor formula
%   Method 2 (Wei ODE45):           Wei's CME via ode45, WITH back-conversion (pump depletion)
%   Method 3 (Brent RK4 / g_Wei):  Brent's RK4, undepleted pump, using Wei's g value
%   Method 4 (Brent RK4 / g_Brent):Brent's RK4, undepleted pump, using Brent's g value
%
% Three figures isolate the two divergences:
%   Fig 1: Wei Analytic vs Wei ODE45       -- internal consistency check
%   Fig 2: Wei ODE45  vs Brent RK4/g_Wei  -- effect of undepleted pump approx (Divergence E)
%   Fig 3: Brent RK4/g_Wei vs /g_Brent    -- effect of overlap factor definition (Divergence B)
%
% DIVERGENCE NOTES mark every point where the two models differ:
%   A -- Refractive indices
%   B -- Overlap factor normalization (critical: ~5x power impact)
%   C -- g formula (algebraically same if overlap is equivalent)
%   D -- Phase singularity workaround in analytic formula
%   E -- Back-conversion in CME (undepleted pump approximation)
%   F -- Sign convention on exp(+-j*dk*z) in SH equation

    close all; clear; clc;

    %% 1. CONSTANTS
    eps0 = 8.854e-12;   % F/m
    c    = 2.9979e8;    % m/s

    %% 2. PARAMETERS — Wei's exact benchmark values
    lam1_m = 0.446e-6;           % pump wavelength (m)
    omega1  = 2*pi*c / lam1_m;   % pump angular frequency (rad/s)

    % DIVERGENCE NOTE A — Refractive indices
    % Wei uses n1 = n2 = 2.04 (same effective index for both pump and SH modes).
    % Brent v2_2 uses n_pump=2.15, n_shg=2.6 (separate COMSOL-derived values).
    % Using Wei's values here for the baseline comparison.
    n1 = 2.04;   % pump effective index
    n2 = 2.04;   % SH effective index (Wei uses same as pump)

    deff = 7.0e-12;   % nonlinear coefficient d33 (m/V)

    % DIVERGENCE NOTE B — Overlap factor normalization  ***CRITICAL***
    % Wei's ol_factor has units of 1/m. It comes from a mode overlap integral whose
    % normalization retains area units:
    %   ol_factor ~ integral(e1^2 * e2 dA) / (norm factors with area units)
    %
    % Brent's overlap_eta is dimensionless (COMSOL integral normalized by each mode
    % independently, so area cancels):
    %   overlap_eta ~ integral(e1^2 * e2 dA) / sqrt(integral(e1^2)^2 * integral(e2^2))
    %
    % Conversion: ol_factor [1/m] = overlap_eta [1] / sqrt(A_eff [m^2])
    %
    % With Brent's geometry:
    %   sqrt(A_eff) = sqrt(300nm * 335nm) = 3.17e-7 m
    %   ol_factor_Brent_equiv = 0.035 / 3.17e-7 = 1.10e5 /m
    %   vs Wei's 4.8e4 /m  =>  ratio ~2.3x  =>  ~5x difference in SHG power
    ol_factor_Wei          = 0.048e6;          % 1/m (Wei's value, TM00+TM04 for W=300nm)
    width_m                = 0.300e-6;         % Brent's waveguide geometry
    height_m               = 0.335e-6;
    A_eff_m2               = width_m * height_m;
    overlap_eta_Brent      = 0.035;            % dimensionless (Brent's COMSOL value)
    ol_factor_Brent_equiv  = overlap_eta_Brent / sqrt(A_eff_m2);  % converted to 1/m

    P_pump_W = 200e-3 * 0.6;   % Wei: 200 mW * 0.6 input coupling = 120 mW

    %% 3. G-COEFFICIENT COMPUTATION
    % DIVERGENCE NOTE C — g formula
    % Both formulas derive from the same nonlinear coupled-mode starting point.
    % They are algebraically identical when ol_factor = overlap_eta / sqrt(A_eff).
    % The numerical difference here comes entirely from Divergence B (overlap normalization).

    % Wei's g: derived from yeta = g^2, using ol_factor in 1/m
    yeta_Wei = 2 * omega1^2 / (eps0 * c^3 * n1^2 * n2) * deff^2 * ol_factor_Wei^2;
    g_Wei_SI = sqrt(yeta_Wei);         % 1/(sqrt(W) * m)
    g_Wei_um = g_Wei_SI * 1e-6;        % 1/(sqrt(W) * um)  -- for Brent's RK4 (z in um)

    % Brent's g_numeric: uses overlap_eta (dimensionless) + A_eff
    g_Brent_SI = (omega1 * deff / c) * sqrt(2 / (eps0 * n1^2 * n2 * c * A_eff_m2)) * overlap_eta_Brent;
    g_Brent_um = g_Brent_SI * 1e-6;   % 1/(sqrt(W) * um)

    fprintf('=== OVERLAP FACTOR COMPARISON ===\n');
    fprintf('  Wei   ol_factor           = %.4e /m\n', ol_factor_Wei);
    fprintf('  Brent ol_factor (equiv)   = %.4e /m  (overlap_eta=%.3f / sqrt(A_eff))\n', ...
            ol_factor_Brent_equiv, overlap_eta_Brent);
    fprintf('  Ratio Brent/Wei           = %.3f\n\n', ol_factor_Brent_equiv / ol_factor_Wei);
    fprintf('=== G-COEFFICIENT COMPARISON ===\n');
    fprintf('  g_Wei   = %.4e 1/(sqrt(W)*m)\n', g_Wei_SI);
    fprintf('  g_Brent = %.4e 1/(sqrt(W)*m)\n', g_Brent_SI);
    fprintf('  Ratio g_Brent/g_Wei = %.3f  =>  SHG power ratio = %.2fx\n\n', ...
            g_Brent_SI/g_Wei_SI, (g_Brent_SI/g_Wei_SI)^2);

    %% 4. SWEEP SETUP
    pump_losses_dBcm = [5, 10, 20, 35];   % dB/cm pump loss values to sweep
    num_losses       = length(pump_losses_dBcm);

    % Length sweep: Wei's range 0.05 mm to 2.0 mm
    L_mm = linspace(0.05, 2.0, 80);
    L_m  = L_mm * 1e-3;   % meters  -- for Wei's ODE45 and analytic formula
    L_um = L_mm * 1e3;    % microns -- for Brent's RK4 engine

    %% 5. ALLOCATE RESULT MATRICES  [length_idx x loss_idx]
    P_Wei_analytic     = zeros(length(L_mm), num_losses);
    P_Wei_ode45        = zeros(length(L_mm), num_losses);
    P_Brent_RK4_gWei   = zeros(length(L_mm), num_losses);
    P_Brent_RK4_gBrent = zeros(length(L_mm), num_losses);

    % DIVERGENCE NOTE D — phase singularity workaround
    % Wei's closed-form PM_factor has a 0/0 singularity when dk=0 AND dalpha=0 simultaneously.
    % Wei avoids it by setting dbeta=1 rad/m instead of zero.
    % At 1 rad/m over 2mm, sin(dbeta*L/2)^2 << 1, so the effect on results is negligible.
    % Brent's RK4 uses dk=0 with no singularity.
    dbeta_Wei = 1;   % rad/m (Wei's singularity workaround — not a physical phase mismatch)

    %% 6. COMPUTE ALL FOUR METHODS
    fprintf('Running simulations...\n');
    for lc = 1:num_losses
        a_pump_dBcm = pump_losses_dBcm(lc);
        a_SH_dBcm   = 10 * a_pump_dBcm;    % SH loss = 10x pump loss (matches Wei and Brent default)

        % --- Convert losses for Wei's code (Np/m) ---
        % Wei multiplies dB/cm by 100 to get dB/m, then divides by 4.343 to get Np/m
        alpha1 = (a_pump_dBcm * 100) / 4.343;   % Np/m
        alpha2 = (a_SH_dBcm   * 100) / 4.343;   % Np/m
        dalpha = alpha2/2 - alpha1;

        % --- Convert losses for Brent's RK4 (Np/um) ---
        % dB/cm / 4.3429 gives Np/cm; * 1e-4 converts to Np/um
        a0_um = (a_pump_dBcm / 4.3429) * 1e-4;  % Np/um
        as_um = (a_SH_dBcm   / 4.3429) * 1e-4;  % Np/um

        % --- Method 1: Wei Analytic ---
        % Closed-form solution to the SHG CME with asymmetric loss and phase mismatch.
        % Reference: Wei's code lines 57-62.
        PM = exp(-alpha2 * L_m) .* ...
             ( (exp(dalpha*L_m) - 1).^2 + 4*exp(dalpha*L_m).*(sin(dbeta_Wei*L_m/2)).^2 ) ...
             ./ ( (dalpha*L_m).^2 + (dbeta_Wei*L_m).^2 );
        P_Wei_analytic(:, lc) = yeta_Wei * L_m.^2 .* PM * P_pump_W^2;

        % --- Method 2: Wei ODE45 (full CME with back-conversion) ---
        % DIVERGENCE NOTE E — back-conversion (pump depletion)
        % Wei's CME for the pump field:
        %   dA1/dz = +j*g*conj(A1)*A2*exp(+j*dk*z) - alpha1/2*A1   <-- back-conversion term present
        % Brent's RK4 (shg_derivs):
        %   dA1/dz = -alpha1/2 * A1                                  <-- undepleted pump, no back-conversion
        %
        % DIVERGENCE NOTE F — sign of phase term in SH equation
        % Wei:   dA2/dz = j*g*A1^2*exp(-j*dk*z) - alpha2/2*A2   (negative sign)
        % Brent: dA3/dz = j*g*A1^2*exp(+j*dk*z) - alpha3/2*A3   (positive sign)
        % At dk=0 (perfect phase match) this makes no difference.
        % At nonzero dk the two models accumulate phase in opposite directions.
        fun_Wei = @(z, x) [ 1i*g_Wei_SI*conj(x(1))*x(2)*exp(+1i*dbeta_Wei*z) - alpha1/2*x(1); ...
                             1i*g_Wei_SI*x(1)^2    *exp(-1i*dbeta_Wei*z)        - alpha2/2*x(2) ];
        x0   = [sqrt(P_pump_W); 0];
        opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-12);
        % Prepend z=0 so ode45 starts from the initial condition
        [~, B] = ode45(fun_Wei, [0, L_m], x0, opts);
        P_Wei_ode45(:, lc) = abs(B(2:end, 2)).^2;   % skip z=0 row

        % --- Methods 3 & 4: Brent RK4 (undepleted pump) ---
        for li = 1:length(L_um)
            [Pg_gW, ~] = rk4_engine(L_um(li), 0.25, g_Wei_um,   a0_um, as_um, as_um, 0, P_pump_W);
            [Pg_gB, ~] = rk4_engine(L_um(li), 0.25, g_Brent_um, a0_um, as_um, as_um, 0, P_pump_W);
            P_Brent_RK4_gWei(li,   lc) = Pg_gW;
            P_Brent_RK4_gBrent(li, lc) = Pg_gB;
        end

        fprintf('  Loss case %d/%d done (pump=%d dB/cm, SH=%d dB/cm)\n', ...
                lc, num_losses, a_pump_dBcm, a_SH_dBcm);
    end

    %% 7. PLOTS
    colors      = lines(num_losses);
    loss_labels = arrayfun(@(x) sprintf('\\alpha_0=%d dB/cm', x), pump_losses_dBcm, 'UniformOutput', false);
    mk_idx      = 1:8:length(L_mm);   % sparse marker positions for overlay visibility

    % --- Figure 1: Wei Internal Consistency (Analytic vs ODE45) ---
    % These two should agree closely. Any gap reveals inconsistency within Wei's own model.
    figure('Color','w','Name','Fig 1: Wei Internal Consistency','Position',[50,600,750,480]);
    hold on;
    h1 = gobjects(num_losses, 1);
    for lc = 1:num_losses
        h1(lc) = plot(L_mm, P_Wei_analytic(:,lc)*1e6, '-', 'Color', colors(lc,:), 'LineWidth', 2.0);
                 plot(L_mm(mk_idx), P_Wei_ode45(mk_idx,lc)*1e6, 'o', 'Color', colors(lc,:), ...
                      'MarkerSize', 5, 'MarkerFaceColor', 'w', 'LineWidth', 1.2);
    end
    set(gca, 'YScale', 'log'); grid on;
    xlabel('Waveguide Length (mm)'); ylabel('SHG Power (\muW)');
    title({'Fig 1: Wei Internal Consistency — Analytic vs ODE45'; ...
           'Solid lines = Wei Analytic,  Open circles = Wei ODE45 (back-conversion included)'});
    legend(h1, loss_labels, 'Location', 'northwest', 'FontSize', 8);

    % --- Figure 2: Back-Conversion Effect ---
    % Both curves use g_Wei and the same losses.
    % Difference is solely due to the undepleted pump approximation (Divergence E).
    % Open circles (Brent RK4) are plotted last so they sit on top of the solid lines (Wei ODE45).
    figure('Color','w','Name','Fig 2: Back-Conversion Effect','Position',[100,100,750,480]);
    hold on;
    h2 = gobjects(num_losses, 1);
    for lc = 1:num_losses
        h2(lc) = plot(L_mm, P_Wei_ode45(:,lc)*1e6, '-', 'Color', colors(lc,:), 'LineWidth', 2.5);
                 plot(L_mm, P_Brent_RK4_gWei(:,lc)*1e6, '--', 'Color', colors(lc,:), 'LineWidth', 1.2);
                 plot(L_mm(mk_idx), P_Brent_RK4_gWei(mk_idx,lc)*1e6, 'o', 'Color', colors(lc,:), ...
                      'MarkerSize', 6, 'MarkerFaceColor', 'w', 'LineWidth', 1.4);
    end
    set(gca, 'YScale', 'log'); grid on;
    xlabel('Waveguide Length (mm)'); ylabel('SHG Power (\muW)');
    title({'Fig 2: Effect of Undepleted Pump Approximation  (same g_{Wei}, same losses)'; ...
           'Solid = Wei ODE45 (back-conversion),  Dashed + Open circles = Brent RK4 (undepleted pump)'});
    legend(h2, loss_labels, 'Location', 'northwest', 'FontSize', 8);

    % --- Figure 3: Overlap Factor (g) Effect ---
    % Both curves use Brent's RK4 with undepleted pump.
    % Difference is solely from the overlap normalization (Divergence B/C).
    figure('Color','w','Name','Fig 3: Overlap Factor Effect','Position',[150,300,750,480]);
    hold on;
    h3 = gobjects(num_losses, 1);
    for lc = 1:num_losses
        h3(lc) = plot(L_mm, P_Brent_RK4_gWei(:,lc)*1e6, '-', 'Color', colors(lc,:), 'LineWidth', 2.5);
                 plot(L_mm, P_Brent_RK4_gBrent(:,lc)*1e6, '--', 'Color', colors(lc,:), 'LineWidth', 1.2);
                 plot(L_mm(mk_idx), P_Brent_RK4_gBrent(mk_idx,lc)*1e6, 'o', 'Color', colors(lc,:), ...
                      'MarkerSize', 6, 'MarkerFaceColor', 'w', 'LineWidth', 1.4);
    end
    set(gca, 'YScale', 'log'); grid on;
    xlabel('Waveguide Length (mm)'); ylabel('SHG Power (\muW)');
    title({sprintf('Fig 3: Overlap Factor Effect  (both Brent RK4, undepleted pump)'); ...
           sprintf('Solid = g_{Wei}  (ol=%.3ge6/m),  Dashed + Open circles = g_{Brent}  (\\eta=%.3f)', ...
                   ol_factor_Wei/1e6, overlap_eta_Brent)});
    legend(h3, loss_labels, 'Location', 'northwest', 'FontSize', 8);

    %% 8. SUMMARY TABLE
    L_ref_mm = 1.0;
    fprintf('\n=== SHG POWER AT L = %.1f mm (uW) ===\n', L_ref_mm);
    fprintf('%-32s', 'Method');
    for lc = 1:num_losses
        fprintf('| a0=%2ddB/cm ', pump_losses_dBcm(lc));
    end
    fprintf('\n%s\n', repmat('-', 32 + num_losses*13, 1));
    method_names = {'Wei Analytic', 'Wei ODE45 (back-conv.)', ...
                    'Brent RK4 (g_Wei)', 'Brent RK4 (g_Brent)'};
    datasets = {P_Wei_analytic, P_Wei_ode45, P_Brent_RK4_gWei, P_Brent_RK4_gBrent};
    for m = 1:4
        fprintf('%-32s', method_names{m});
        for lc = 1:num_losses
            val = interp1(L_mm, datasets{m}(:,lc), L_ref_mm) * 1e6;
            fprintf('| %9.3f uW ', val);
        end
        fprintf('\n');
    end
    fprintf('%s\n\n', repmat('-', 32 + num_losses*13, 1));
    fprintf('NOTE: At low efficiency (current regime), Methods 2 and 3 should be nearly identical.\n');
    fprintf('      Any gap between them quantifies the back-conversion error (Divergence E).\n');
    fprintf('      The gap between Methods 3 and 4 quantifies the overlap factor issue (Divergence B).\n');
end


%% =========================================================================
%  RK4 ENGINE  (from Brent SHG_Design_Suite_v2_2)
%  State vector Y = [Re(A1), Im(A1), Re(A3), Im(A3), Ps]
%  A1 = pump field amplitude (sqrt(W)), A3 = SH field amplitude, Ps = scattered SH power (W)
%  All z-units in um; g in 1/(sqrt(W)*um); a0, a3t, a3s in Np/um
% =========================================================================
function [Pg_end, Ps_end] = rk4_engine(L, dz, g, a0, a3t, a3s, dk, Pp)
    z_vec = 0:dz:L;
    Y = [sqrt(Pp); 0; 0; 0; 0];
    for i = 1:(length(z_vec)-1)
        z  = z_vec(i);
        k1 = shg_derivs(z,       Y,          a0, a3t, a3s, g, dk);
        k2 = shg_derivs(z+dz/2, Y+k1*dz/2,  a0, a3t, a3s, g, dk);
        k3 = shg_derivs(z+dz/2, Y+k2*dz/2,  a0, a3t, a3s, g, dk);
        k4 = shg_derivs(z+dz,   Y+k3*dz,    a0, a3t, a3s, g, dk);
        Y  = Y + (dz/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
    Pg_end = Y(3)^2 + Y(4)^2;
    Ps_end = Y(5);
end

function dY = shg_derivs(z, Y, a0, a3t, a3s, g, dk)
    A1 = Y(1) + 1i*Y(2);
    A3 = Y(3) + 1i*Y(4);
    dA1 = -0.5*a0*A1;                                      % undepleted pump (Divergence E)
    dA3 = -0.5*a3t*A3 + 1i*g*(A1^2)*exp(+1i*dk*z);        % +sign on dk (Divergence F)
    dPs = a3s * (abs(A3)^2);
    dY  = [real(dA1); imag(dA1); real(dA3); imag(dA3); dPs];
end

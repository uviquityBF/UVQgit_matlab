%% MATLAB Simulation: Short-Delay Self-Heterodyne Spectrum (Jin et al., 2022)
% Reproduces Figure 2(a), 2(b), and 2(c) from:
% Duo Jin et al., "Narrow laser-linewidth measurement using short delay
% self-heterodyne interferometry," Optics Express 30(17), 30600 (2022).

clear; clc; close all;

%% ========================================================================
%  PHYSICAL CONSTANTS & CORE FUNCTIONS
%  ========================================================================
c = 3e8;            % Speed of light in vacuum (m/s)
n_default = 1.46;   % Refractive index of fiber core

% Function handle for S2(f) [Equation 3 in Jin et al.]
S2_func = @(f, f1, df, tau_d) ...
    1 - exp(-2*pi*tau_d*df) .* ( cos(2*pi*tau_d*(f - f1)) + ...
    df .* sinc_term(f, f1, tau_d) );

% Function handle for S1(f) [Equation 2 in Jin et al.]
S1_func = @(f, f1, df, P0) ...
    (P0^2 / (4*pi)) .* (df ./ (df^2 + (f - f1).^2));

% Function handle for Total PSD S(f) [Equation 1, ignoring delta S3 off-center]
S_func = @(f, f1, df, tau_d, P0) S1_func(f, f1, df, P0) .* S2_func(f, f1, df, tau_d);

%% ========================================================================
%  PART 1: REPRODUCE FIGURE 2 FROM JIN ET AL. (2022)
%  ========================================================================
f1_paper = 100e6;      % AOM Modulation Frequency = 100 MHz
P0_paper = 1;          % Optical Power scaling factor
n = n_default;

% Frequency vector spanning 99.5 MHz to 100.5 MHz (1 MHz total span)
f_paper = linewidth_span(f1_paper, 1e6, 10001);

% --- Figure 2(a): Fixed linewidth (1 kHz), varying fiber lengths ---
df_fig2a = 1e3; % 1 kHz
L_vec_2a = [1e3, 3e3, 10e3, 50e3, 500e3]; % Fiber lengths: 1km, 3km, 10km, 50km, 500km
legend_2a = {'L = 1 km', 'L = 3 km', 'L = 10 km', 'L = 50 km', 'L = 500 km'};

figure('Name', 'Figure 2(a) - Varying Fiber Lengths', 'Color', 'w');
hold on;
for i = 1:length(L_vec_2a)
    L_val = L_vec_2a(i);
    tau_d = (n * L_val) / c;
    S_vals = S_func(f_paper, f1_paper, df_fig2a, tau_d, P0_paper);

    % Convert to relative dB scale
    S_dB = 10 * log10(S_vals);
    plot(f_paper / 1e6, S_dB, 'LineWidth', 1.5);
end
hold off;
grid on;
title('Fig 2(a): PSD for Fixed \Deltaf = 1 kHz, Varying Delay Fiber Length L');
xlabel('Frequency (MHz)');
ylabel('Power Spectrum (a.u. / dB)');
legend(legend_2a, 'Location', 'northeast');

% --- Figure 2(b): Fixed fiber length (3 km), varying linewidths ---
L_fig2b = 3e3; % 3 km
tau_d_2b = (n * L_fig2b) / c;
df_vec_2b = [0.1e3, 1e3, 5e3, 10e3, 50e3]; % 0.1kHz, 1kHz, 5kHz, 10kHz, 50kHz
legend_2b = {'\Deltaf = 0.1 kHz', '\Deltaf = 1 kHz', '\Deltaf = 5 kHz', ...
             '\Deltaf = 10 kHz', '\Deltaf = 50 kHz'};

figure('Name', 'Figure 2(b) - Varying Linewidths', 'Color', 'w');
hold on;
for i = 1:length(df_vec_2b)
    df_val = df_vec_2b(i);
    S_vals = S_func(f_paper, f1_paper, df_val, tau_d_2b, P0_paper);

    S_dB = 10 * log10(S_vals);
    plot(f_paper / 1e6, S_dB, 'LineWidth', 1.5);
end
hold off;
grid on;
title('Fig 2(b): PSD for Fixed L = 3 km, Varying Laser Linewidth \Deltaf');
xlabel('Frequency (MHz)');
ylabel('Power Spectrum (a.u. / dB)');
legend(legend_2b, 'Location', 'northeast');

% --- Figure 2(c): Function S2(f) for different delay fiber lengths ---
figure('Name', 'Figure 2(c) - S2 Function Envelope', 'Color', 'w');
for i = 1:length(L_vec_2a)
    L_val = L_vec_2a(i);
    tau_d = (n * L_val) / c;
    S2_vals = S2_func(f_paper, f1_paper, df_fig2a, tau_d);

    subplot(length(L_vec_2a), 1, i);
    plot(f_paper / 1e6, S2_vals, 'LineWidth', 1.2, 'Color', [0 0.447 0.741]);
    grid on;
    ylabel('S_2 Amplitude');
    title(sprintf('L = %d km', L_val/1e3));
    ylim([-0.2 2.2]);
end
xlabel('Frequency (MHz)');

%% ========================================================================
%  PART 2: PREDICTED TRACE(S) FOR YOUR EXPERIMENT
%  ========================================================================
% my_df and/or my_L may each be a scalar or a vector. If one is a vector
% and the other a scalar, the scalar is held fixed while the vector is
% swept. If both are vectors, they must be the same length and are paired
% element-wise. Every case is overlaid on one plot.

% --- YOUR EXPERIMENTAL INPUT PARAMETERS ---
my_df   = [1e3, 1e4, 1e5, 1e6, 10e6];  % Candidate laser linewidths (Hz)
my_L    = 20;             % Fiber Delay Length (m)
my_f1   = 80e6;          % AOM Modulation Frequency: 80 MHz
my_n    = 1.46;          % PM405 Fiber Refractive Index
my_P0   = 1;             % Normalized Optical Power

nCases = max(numel(my_df), numel(my_L));
if ~(numel(my_df) == 1 || numel(my_df) == nCases)
    error('my_df must be scalar or the same length as my_L.');
end
if ~(numel(my_L) == 1 || numel(my_L) == nCases)
    error('my_L must be scalar or the same length as my_df.');
end

% Size the frequency domain so it reaches at least the first two extrema
% (Eq. 7, m=0 and m=1) for every case, regardless of fiber length.
margin = 1.3;
df_m1_allcases = zeros(1, nCases);
for k = 1:nCases
    Lk = my_L(min(k, numel(my_L)));
    df_m1_allcases(k) = (3*c) / (2 * my_n * Lk);
end
half_span = max(max(df_m1_allcases) * margin, 5 * max(my_df));
f_exp = linewidth_span(my_f1, 2 * half_span, 20001);

fprintf('\n================ YOUR EXPERIMENTAL PREDICTIONS ================\n');
figure('Name', 'Your Experiment - Predicted Scope Trace', 'Color', 'w');
hold on;
colors = lines(nCases);

for k = 1:nCases
    dfk = my_df(min(k, numel(my_df)));
    Lk  = my_L(min(k, numel(my_L)));
    tau_dk = (my_n * Lk) / c;

    % Theoretical offsets for first valley (m=0) and first peak (m=1),
    % reported for comparison only -- see Part 3 for when these hold.
    df_m0 = c / (my_n * Lk);
    df_m1 = (3*c) / (2 * my_n * Lk);

    S_exp = S_func(f_exp, my_f1, dfk, tau_dk, my_P0);
    S_exp_dB = 10 * log10(S_exp);

    % Locate the true first minimum and first maximum numerically (valid
    % even when Delta_f is comparable to the fringe spacing, unlike Eq. 7).
    upperMask = f_exp > my_f1;
    f_upper = f_exp(upperMask);
    S_upper = S_exp_dB(upperMask);
    idxMin = find(islocalmin(S_upper), 1, 'first');
    idxMaxAll = find(islocalmax(S_upper));

    fprintf('--- Case %d: Delta_f = %.4g Hz, L = %.4g m ---\n', k, dfk, Lk);
    fprintf('Eq. 7 estimate: valley at +%.3f MHz, peak at +%.3f MHz\n', ...
            df_m0/1e6, df_m1/1e6);

    if isempty(idxMin)
        warning('Case %d: no resolvable local minimum -- fringes are washed out (pure Lorentzian).', k);
        plot(f_exp/1e6, S_exp_dB, 'Color', colors(k,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('\\Deltaf=%.3g Hz, L=%.3g m', dfk, Lk));
        continue;
    end
    idxMaxAfter = idxMaxAll(idxMaxAll > idxMin);
    if isempty(idxMaxAfter)
        warning('Case %d: no resolvable local maximum after the first valley.', k);
        idxMax = [];
    else
        idxMax = idxMaxAfter(1);
    end

    fL = f_upper(idxMin); SL_dB = S_upper(idxMin);
    plot(f_exp/1e6, S_exp_dB, 'Color', colors(k,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('\\Deltaf=%.3g Hz, L=%.3g m', dfk, Lk));
    plot(fL/1e6, SL_dB, 'v', 'Color', colors(k,:), 'MarkerFaceColor', colors(k,:), ...
         'HandleVisibility', 'off');

    if ~isempty(idxMax)
        fH = f_upper(idxMax); SH_dB = S_upper(idxMax);
        plot(fH/1e6, SH_dB, '^', 'Color', colors(k,:), 'MarkerFaceColor', colors(k,:), ...
             'HandleVisibility', 'off');
        deltaS = SH_dB - SL_dB;
        fprintf('Numeric result: valley at +%.3f MHz (%.2f dB), peak at +%.3f MHz (%.2f dB)\n', ...
                (fL-my_f1)/1e6, SL_dB, (fH-my_f1)/1e6, SH_dB);
        fprintf('Numeric contrast (Delta S): %.2f dB\n', deltaS);
    end
end
hold off;
grid on;
xlabel('Frequency (MHz)');
ylabel('Power Spectral Density (dB)');
title(sprintf('Predicted Short-DSH Spectrum (f_1 = %d MHz) -- triangle-down = valley, triangle-up = peak', my_f1/1e6));
legend('Location', 'eastoutside');
fprintf('=================================================================\n');

%% ========================================================================
%  PART 3: REGIME MAP -- FRINGE VISIBILITY vs FIBER LENGTH & LINEWIDTH
%  ========================================================================
% The coherence-envelope oscillations are damped by the exp(-2*pi*tau_d*Delta_f)
% term in S2_func, which depends only on the product x = tau_d*Delta_f =
% n*L*Delta_f/c. Small x -> strong fringes; x >> 1 -> fringes wash out and
% S(f) reduces to a pure Lorentzian. This map shows that visibility,
% V = exp(-2*pi*x), across a broad range of fiber length and linewidth, so
% you can see directly whether shortening the fiber recovers the fringes
% for a given linewidth.

L_grid  = logspace(-1, 2, 200);   % 0.1 m to 100 m
df_grid = logspace(4, 8, 200);    % 10 kHz to 100 MHz
[L_mesh, df_mesh] = meshgrid(L_grid, df_grid);

x_mesh = (my_n .* L_mesh .* df_mesh) / c;
V_mesh = exp(-2*pi*x_mesh);

figure('Name', 'Regime Map - Fringe Visibility', 'Color', 'w');
contourf(L_mesh, df_mesh/1e6, V_mesh, 20, 'LineStyle', 'none');
set(gca, 'XScale', 'log', 'YScale', 'log');
colormap(flipud(parula));
cb = colorbar;
cb.Label.String = 'Fringe visibility  V = exp(-2\pi\tau_d\Deltaf)';
hold on;

% Reference contour marking a practically "usable" visibility threshold
[C, hThresh] = contour(L_mesh, df_mesh/1e6, V_mesh, [0.3 0.3], 'k--', 'LineWidth', 1.5);
clabel(C, hThresh, 'FontSize', 9);

% Recommended-length curve for a target visibility of x = 0.2
L_opt = 0.2 * c ./ (my_n .* df_grid);
hOpt = plot(L_opt, df_grid/1e6, 'w-', 'LineWidth', 1.5);

% This script's simulated Part 2 cases, for reference
caseL  = arrayfun(@(k) my_L(min(k, numel(my_L))), 1:nCases);
caseDf = arrayfun(@(k) my_df(min(k, numel(my_df))), 1:nCases);
hCases = plot(caseL, caseDf/1e6, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);

xlabel('Delay Fiber Length, L (m)');
ylabel('Laser Linewidth, \Deltaf (MHz)');
title('Regime Map: Coherence-Envelope Fringe Visibility');
legend([hThresh, hOpt, hCases], {'V = 0.3 threshold', 'L for target x=0.2', 'This script''s cases'}, ...
       'Location', 'northoutside');

%% ========================================================================
%  HELPER FUNCTIONS FOR GRACEFUL NUMERICAL EVALUATION
%  ========================================================================

% Evaluates sin(x)/x without division-by-zero at f = f1
function term = sinc_term(f, f1, tau_d)
    x = f - f1;
    term = zeros(size(x));

    % For x != 0: sin(2*pi*tau_d*x) / x
    idx_non_zero = (x ~= 0);
    term(idx_non_zero) = sin(2*pi*tau_d * x(idx_non_zero)) ./ x(idx_non_zero);

    % Limit as x -> 0 is 2*pi*tau_d
    term(~idx_non_zero) = 2*pi*tau_d;
end

% Generates a symmetric frequency vector centered at f1
function f_vec = linewidth_span(f1, span, num_points)
    f_vec = linspace(f1 - span/2, f1 + span/2, num_points);
end

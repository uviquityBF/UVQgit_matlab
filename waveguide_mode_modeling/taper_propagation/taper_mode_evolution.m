%% Mode Evolution Down a Taper (eigenmode expansion)
% What happens to the horn modes as the waveguide narrows? Which survive to the
% narrow end, which are radiated (cut off), and which get converted into the
% fundamental or into other modes?
%
% Method: the taper is a staircase of sections centred on the exported widths.
% At each junction the modes of one width are projected onto the modes of the
% next (overlap matrix, E-field transverse inner product); within a section each
% mode accumulates phase exp(i*2*pi*neff*L/lambda) (and attenuation if its COMSOL
% neff has an imaginary part). Forward propagation only, no reflections. Modes
% with Re(neff) below the substrate/cladding index are treated as radiated.
%
% Inputs
%   - a directory of combined COMSOL mode exports, one file per width, with
%     "w=<number>nm" in the file name (modeData_N=20__h=350nm__w=2500nm.csv).
%     All files must be on the same regular x,y grid.
%   - the input light, either the Zemax spot (decomposed into horn modes exactly as in
%     modal_decomposition/modal_power_distribution.m), the fundamental only, or equal
%     power in all horn modes (cfg.input_source).
% Taper shape and length are free parameters: the overlap matrices are computed once and
% reused for every profile / length. Compatible with MATLAB R2018b.

clearvars -except cfg_override; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'helpers', 'waveguide'));

%% ---------------- USER SETTINGS ----------------
cfg = modal_default_cfg();
cfg.grid_n = 501;                    % grid for the spot decomposition (input_source = 'spot')

cfg.mode_dir       = '';             % '' = ask; directory of per-width mode files
cfg.n_substrate    = 1.8;            % sapphire
cfg.n_clad_top     = 1.48;           % SiO2
cfg.drop_leaky     = true;           % true: Re(neff) < max(n_sub, n_clad) => radiated (removed)
                                     % false: keep them, using Im(neff) as attenuation
cfg.lambda_mode_um = 0.45;           % wavelength of the COMSOL modes (for the phases)
cfg.crop_x_um      = [-5 5];         % crop window for the overlaps (um)
cfg.crop_y_um      = [-3 3];
cfg.w_start_um     = [2.5];          % taper starting width;   [] = widest exported width
cfg.w_end_um       = [0.7];          % taper ending width;   [] = narrowest exported width

cfg.input_source   = 'equal';         % 'spot' | 'fundamental' | 'equal'
cfg.lengths_um     = [25 50 100 250 500 1000 2000];   % taper lengths for the length sweep
cfg.L_detail_um    = 350;            % length used for the detailed plots
cfg.profile_detail = 1;              % which row of cfg.profiles is used for the detailed plots
% name, f(s): remaining fraction of the width change, f(0)=1, f(1)=0, decreasing
cfg.profiles = { ...
    'linear',                   @(s) 1 - s;
    'fast narrowing first',     @(s) (1 - s).^2;
    'slow narrowing first',     @(s) 1 - s.^2;
    'S-curve (cosine)',         @(s) 0.5 * (1 + cos(pi * s)) };

cfg.check_stride   = 2;              % convergence check: also run with every n-th width ([] = skip)
cfg.mc_trials      = 0;              % >0: Monte-Carlo test of sensitivity to neff rounding (header has 4 digits)
cfg.neff_jitter    = 5e-5;           % +/- uniform jitter applied to every neff in that test

if exist('cfg_override', 'var')      % optional: struct of settings set before running (batch use)
    fn = fieldnames(cfg_override);
    for i = 1:numel(fn), cfg.(fn{i}) = cfg_override.(fn{i}); end
end
um = char(181);
n_cut = max(cfg.n_substrate, cfg.n_clad_top);

%% ---------------- 1. LOAD MODE SETS (one per width) ----------------
fprintf('--- Step 1: mode sets ---\n');
if isempty(cfg.mode_dir)
    start = pwd;  if exist(cfg.default_dir_modes, 'dir'), start = cfg.default_dir_modes; end
    cfg.mode_dir = uigetdir(start, 'Select directory of per-width COMSOL mode files');
    if isequal(cfg.mode_dir, 0), disp('Canceled.'); return; end
end
opts = struct('n_cut', n_cut, 'drop_leaky', cfg.drop_leaky, ...
              'crop_x_um', cfg.crop_x_um, 'crop_y_um', cfg.crop_y_um);
[W, ~] = load_mode_width_series(cfg.mode_dir, opts);

if ~isempty(cfg.w_start_um), W = W([W.width_um] <= cfg.w_start_um + 1e-9); end
if ~isempty(cfg.w_end_um),   W = W([W.width_um] >= cfg.w_end_um - 1e-9);   end
Nw = numel(W);
if Nw < 2, error('Need at least two widths in %s (found %d).', cfg.mode_dir, Nw); end

fprintf('\n%-9s %-6s %-8s %-6s %-10s %-10s\n', 'width(um)', 'modes', 'dropped', 'kept', 'neff(max)', 'centroid(x,y)');
for i = 1:Nw
    fprintf('%-9.3f %-6d %-8d %-6d %-10.4f (%.2f, %.2f)\n', W(i).width_um, W(i).n_all, W(i).n_dropped, ...
            numel(W(i).neff), W(i).neff(1), W(i).cx_um, W(i).cy_um);
end
if max(abs([W.cx_um] - W(1).cx_um)) > 0.3 || max(abs([W.cy_um] - W(1).cy_um)) > 0.3
    fprintf(2, 'WARNING: fundamental-mode centroids differ between widths by > 0.3 um; are the waveguides centred consistently in COMSOL?\n');
end
if any([W.n_dropped] > 0)
    fprintf('Modes below n = %.3f are treated as radiated (cfg.drop_leaky).\n', n_cut);
end

%% ---------------- 2. INPUT AMPLITUDES AT THE WIDEST SECTION ----------------
fprintf('\n--- Step 2: input amplitudes (%s) ---\n', cfg.input_source);
k1 = numel(W(1).neff);
tags = arrayfun(@(k) sprintf('m%02d', k), (1:k1)', 'UniformOutput', false);
switch lower(cfg.input_source)
    case 'spot'
        cfgS = cfg;  cfgS.mode_file = fullfile(cfg.mode_dir, W(1).file);
        P = load_modal_problem(cfgS);
        if isempty(P), return; end
        if P.nM ~= W(1).n_all, error('Mode count mismatch between spot decomposition and width series.'); end
        shift = cfg.beam_offset_um;  if cfg.align_centroids, shift = shift + P.align_shift; end
        E_in = modal_input_field(P, shift, cfg.defocus_um, cfg.extra_phase_fn);
        nx = numel(P.intX);  ny = numel(P.intY);
        wx = ones(1, nx) * (P.intX(2) - P.intX(1));  wx([1 end]) = wx([1 end]) / 2;
        wy = ones(ny, 1) * (P.intY(2) - P.intY(1));  wy([1 end]) = wy([1 end]) / 2;
        Wt = wy * wx;
        Mm = reshape(P.E_modes, [], P.nM);  Mo = reshape(P.E_other, [], P.nM);
        Pm = (abs(Mm).^2 + abs(Mo).^2)' * Wt(:);
        a0_all = (Mm' * (Wt(:) .* E_in(:))) ./ sqrt(Pm * P.P_in_total);
        a0 = a0_all(W(1).keep_idx);
        tags = P.modeTag(W(1).keep_idx);
        fprintf('Power captured by the kept horn modes: %.2f %% of the incident spot\n', 100 * sum(abs(a0).^2));
    case 'fundamental'
        a0 = zeros(k1, 1);  a0(1) = 1;
    case 'equal'
        a0 = ones(k1, 1) / sqrt(k1);
    otherwise
        error('cfg.input_source must be ''spot'', ''fundamental'' or ''equal''.');
end

%% ---------------- 3. OVERLAP MATRICES (computed once) ----------------
fprintf('\n--- Step 3: junction overlaps ---\n');
chain = make_chain(W, 1:Nw);
if ~isempty(cfg.check_stride) && cfg.check_stride > 1 && Nw >= 3
    idxC = unique([1:cfg.check_stride:Nw, Nw]);
    chainC = make_chain(W, idxC);
else
    chainC = [];
end
% how much of each mode's power survives each junction (sanity / adiabaticity info)
for j = 1:Nw - 1
    cap = sum(abs(chain.C{j}).^2, 1);        % per OLD mode: power projected onto the guided NEW modes
    fprintf('  %.3f -> %.3f um : fundamental keeps %.2f %% as fundamental, %.2f %% in guided modes (modes 1..%d: min %.2f %%)\n', ...
        W(j).width_um, W(j + 1).width_um, 100 * abs(chain.C{j}(1, 1))^2, 100 * cap(1), ...
        min(4, numel(cap)), 100 * min(cap(1:min(4, numel(cap)))));
end

%% ---------------- 4. RUN: all profiles x lengths ----------------
nP = size(cfg.profiles, 1);  nL = numel(cfg.lengths_um);
P_fund = zeros(nP, nL);  P_ho = zeros(nP, nL);  P_tot = zeros(nP, nL);
for p = 1:nP
    for l = 1:nL
        [~, a_out] = run_taper(chain, a0, cfg.profiles{p, 2}, cfg.lengths_um(l), cfg.lambda_mode_um);
        pw = abs(a_out).^2;
        P_fund(p, l) = pw(1);  P_ho(p, l) = sum(pw(2:end));  P_tot(p, l) = sum(pw);
    end
end
P0_fund = abs(a0(1))^2;  P0_ho = sum(abs(a0(2:end)).^2);
if max(P_tot(:)) > 1.02 * sum(abs(a0).^2)
    fprintf(2, 'WARNING: transmitted power exceeds the input in some cases (E-only overlap is not power-conserving here).\n');
end

fprintf('\nPower at the narrow end, as %% of incident spot (profile: %s)\n', cfg.profiles{cfg.profile_detail, 1});
fprintf('%-10s %-12s %-14s %-12s\n', 'L (um)', 'fundamental', 'higher-order', 'total guided');
fprintf('%-10s %-12.2f %-14.2f %-12.2f   <- at the horn (captured)\n', '0', 100 * P0_fund, 100 * P0_ho, 100 * (P0_fund + P0_ho));
for l = 1:nL
    fprintf('%-10g %-12.2f %-14.2f %-12.2f\n', cfg.lengths_um(l), 100 * P_fund(cfg.profile_detail, l), ...
            100 * P_ho(cfg.profile_detail, l), 100 * P_tot(cfg.profile_detail, l));
end

%% ---------------- 5. DETAIL CASE ----------------
pf = cfg.profiles{cfg.profile_detail, 2};  Ld = cfg.L_detail_um;
[a_sec, a_out_d, z_edges] = run_taper(chain, a0, pf, Ld, cfg.lambda_mode_um);
k_end = numel(W(end).neff);
[~, T] = run_taper(chain, eye(k1), pf, Ld, cfg.lambda_mode_um);      % [k_end x k1] amplitude transmission
Tpow = abs(T).^2;
pw = abs(a_out_d).^2;
fprintf('\nDetail case: %s taper, L = %g %sm. Output mode powers (%% of incident spot):\n', cfg.profiles{cfg.profile_detail, 1}, Ld, um);
for k = 1:min(k_end, 8), fprintf('   narrow-end mode %02d (neff %.4f): %.3f %%\n', k, W(end).neff(k), 100 * pw(k)); end

if ~isempty(chainC)
    [~, aC] = run_taper(chainC, a0, pf, Ld, cfg.lambda_mode_um);
    fprintf('Convergence check (every %d-th width, %d sections vs %d): fundamental %.2f %% vs %.2f %%, total %.2f %% vs %.2f %%\n', ...
        cfg.check_stride, numel(chainC.widths), Nw, 100 * abs(aC(1))^2, 100 * pw(1), 100 * sum(abs(aC).^2), 100 * sum(pw));
end

if cfg.mc_trials > 0
    rng(1);  res = zeros(cfg.mc_trials, 3);
    for t = 1:cfg.mc_trials
        cj = chain;
        for j = 1:Nw, cj.neff_re{j} = chain.neff_re{j} + cfg.neff_jitter * (2 * rand(size(chain.neff_re{j})) - 1); end
        [~, aj] = run_taper(cj, a0, pf, Ld, cfg.lambda_mode_um);
        pj = abs(aj).^2;  res(t, :) = [pj(1), sum(pj(2:end)), sum(pj)];
    end
    fprintf('Sensitivity to neff rounding (+/-%.0e, %d trials): fundamental %.2f +/- %.2f %%, higher-order %.2f +/- %.2f %%, total %.2f +/- %.2f %%\n', ...
        cfg.neff_jitter, cfg.mc_trials, 100 * mean(res(:, 1)), 100 * std(res(:, 1)), 100 * mean(res(:, 2)), ...
        100 * std(res(:, 2)), 100 * mean(res(:, 3)), 100 * std(res(:, 3)));
end

%% ---------------- 6. PLOTS ----------------
col = lines(max(nP, 6));

% Fig 1: neff vs width
figure('Color', 'w', 'Position', [40 40 760 520], 'Name', 'neff vs width');  hold on;  box on;  grid on;
for i = 1:Nw
    plot(W(i).width_um * ones(size(W(i).neff_all)), W(i).neff_all, 'x', 'Color', [0.7 0.7 0.7]);
    plot(W(i).width_um * ones(size(W(i).neff)), W(i).neff, 'o', 'Color', [0.1 0.4 0.8], 'MarkerFaceColor', [0.1 0.4 0.8], 'MarkerSize', 4);
end
plot(xlim, [cfg.n_substrate cfg.n_substrate], 'r--', 'LineWidth', 1.2);
plot(xlim, [cfg.n_clad_top cfg.n_clad_top], 'm--', 'LineWidth', 1.2);
text(max([W.width_um]), cfg.n_substrate, ' substrate', 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
xlabel(['Waveguide width (' um 'm)']);  ylabel('Re(n_{eff})');
title('Mode effective indices vs width  (blue = kept, gray x = below cutoff)');
set(gca, 'XDir', 'reverse');

% Fig 2: transmission matrix (detail case)
figure('Color', 'w', 'Position', [60 60 1100 520], 'Name', 'Transmission matrix');
subplot(1, 2, 1);
imagesc(100 * Tpow);  colorbar;  colormap(gca, parula);
set(gca, 'XTick', 1:k1, 'XTickLabel', tags, 'XTickLabelRotation', 60, 'FontSize', 8);
set(gca, 'YTick', 1:k_end);  ylabel(sprintf('narrow-end mode (width %.2f %sm), 1 = fundamental', W(end).width_um, um));
xlabel(sprintf('horn mode (width %.2f %sm)', W(1).width_um, um));
title(sprintf('Power transfer (%%) | %s, L = %g %sm', cfg.profiles{cfg.profile_detail, 1}, Ld, um));
subplot(1, 2, 2);
bar(1:k1, 100 * sum(Tpow, 1), 0.7, 'FaceColor', [0.15 0.45 0.75], 'EdgeColor', 'none');
set(gca, 'XTick', 1:k1, 'XTickLabel', tags, 'XTickLabelRotation', 60, 'FontSize', 8);
ylabel('Power still guided at narrow end (%)');  grid on;  box on;  ylim([0 105]);
title('Survival of each horn mode');

% Fig 3: along the taper (detail case)
n = numel(a_sec);  kmax = max(cellfun(@(a) size(a, 1), a_sec));
Pz = nan(kmax, n);
for j = 1:n, Pz(1:size(a_sec{j}, 1), j) = abs(a_sec{j}(:, 1)).^2; end
Ptot_z = sum(Pz, 1, 'omitnan');
[~, ordm] = sort(max(Pz, [], 2), 'descend');  show = sort(ordm(1:min(6, kmax)));
figure('Color', 'w', 'Position', [80 80 900 650], 'Name', 'Along the taper');
subplot(2, 1, 1);  hold on;  box on;  grid on;
s = linspace(0, 1, 400);
for p = 1:nP
    lw = 1.2;  if p == cfg.profile_detail, lw = 2.5; end
    plot(Ld * s, W(end).width_um + (W(1).width_um - W(end).width_um) * cfg.profiles{p, 2}(s), 'LineWidth', lw, ...
         'Color', col(p, :), 'DisplayName', cfg.profiles{p, 1});
end
plot((z_edges(1:end-1) + z_edges(2:end)) / 2, [W(1:Nw).width_um], 'k.', 'MarkerSize', 9, 'HandleVisibility', 'off');
ylabel(['Width (' um 'm)']);  legend('Location', 'best', 'FontSize', 8);  title('Taper profiles (dots = centre and width of each section, detailed profile)');
subplot(2, 1, 2);  hold on;  box on;  grid on;
zz = [0 z_edges(2:end)];
for q = 1:numel(show)
    plot(zz, 100 * [abs(a0(min(show(q), k1)))^2 * (show(q) <= k1), Pz(show(q), :)], '-o', 'LineWidth', 1.3, 'MarkerSize', 3, ...
         'DisplayName', sprintf('local mode %02d', show(q)));
end
plot(zz, 100 * [sum(abs(a0).^2), Ptot_z], 'k--', 'LineWidth', 1.5, 'DisplayName', 'total guided');
xlabel(['z along taper (' um 'm)']);  ylabel('Power (% of incident spot)');
title('Power in the local modes (mode index = rank by neff in that section)');  legend('Location', 'best', 'FontSize', 8);

% Fig 4: length sweep
figure('Color', 'w', 'Position', [100 100 1500 420], 'Name', 'Length sweep');
dat = {P_fund, P_ho, P_tot};  ref = [P0_fund, P0_ho, P0_fund + P0_ho];
ttl = {'Fundamental at narrow end', 'Higher-order modes at narrow end', 'Total still guided'};
for q = 1:3
    subplot(1, 3, q);  hold on;  box on;  grid on;
    for p = 1:nP
        plot(cfg.lengths_um, 100 * dat{q}(p, :), '-o', 'LineWidth', 1.5, 'Color', col(p, :), 'DisplayName', cfg.profiles{p, 1});
    end
    plot(xlim, 100 * [ref(q) ref(q)], 'k:', 'LineWidth', 1.2, 'DisplayName', 'at horn (captured)');
    set(gca, 'XScale', 'log');  xlabel(['Taper length (' um 'm)']);  ylabel('% of incident spot');  title(ttl{q});
    if q == 1, legend('Location', 'best', 'FontSize', 8); end
end
sgtitle(sprintf('%s  |  widths %.2f to %.2f %sm (%d sections)  |  input: %s', ...
        cfg.mode_dir, W(1).width_um, W(end).width_um, um, Nw, cfg.input_source), 'Interpreter', 'none', 'FontSize', 8);

if ~isempty(cfg.save_figs_dir)
    figs = findall(0, 'Type', 'figure');
    for fg = figs'
        saveas(fg, fullfile(cfg.save_figs_dir, [regexprep(get(fg, 'Name'), '\W+', '_') '.png']));
    end
end

%% ---------------- local functions ----------------
function chain = make_chain(W, idx)
% Junction matrices between the chosen widths (indices into W, widest first).
    n = numel(idx);
    chain.widths = [W(idx).width_um];
    chain.neff_re = cell(n, 1);  chain.neff_im = cell(n, 1);  chain.C = cell(n - 1, 1);
    for j = 1:n
        chain.neff_re{j} = W(idx(j)).neff;  chain.neff_im{j} = W(idx(j)).neff_imag;
    end
    for j = 1:n - 1
        chain.C{j} = double(W(idx(j + 1)).A' * W(idx(j)).A);     % [k_new x k_old]
    end
end

function [a_sec, a_out, z_edges] = run_taper(chain, a0, profile_fn, L_um, lambda_um)
    [L_sec, z_edges] = taper_section_lengths(chain.widths, profile_fn, L_um);
    [a_sec, a_out] = taper_propagate(a0, chain.C, chain.neff_re, chain.neff_im, L_sec, lambda_um);
end

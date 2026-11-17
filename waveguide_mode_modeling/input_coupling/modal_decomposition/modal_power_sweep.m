%% Modal Power Sweep: beam offset and defocus
% How does the distribution of input power over the waveguide modes change as
% the beam is (a) moved laterally (x offset, y offset, or both) and
% (b) defocused (the POP field is propagated by dz before the overlap)?
%
% Inputs are the same as modal_power_distribution.m (Zemax POP irradiance
% [+ phase], and a combined COMSOL mode .csv for one geometry).
%
% Offsets are relative to the aligned position (input centroid on the
% fundamental-mode centroid) unless cfg.align_centroids = false.
% Positive dz = propagate further along the beam direction (needs phase
% information to be meaningful; flip cfg.phase_sign if Zemax's sign convention
% is opposite). Only a few um of defocus are trustworthy because the POP window
% is small (the part of the beam outside it is dropped).
% Compatible with MATLAB R2018b.

clearvars -except cfg_override sweep_override; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'helpers', 'waveguide'));

%% ---------------- USER SETTINGS ----------------
cfg = modal_default_cfg();
cfg.spot_source = 'cps532';

cfg.grid_n = 501;                    % coarser than the single-run script: speed (sweeps loop a lot)

sweep.dx_um     = -3:0.25:3;         % 1-D sweep of x offset  ([] to skip)
sweep.dy_um     = -1:0.1:1;          % 1-D sweep of y offset  ([] to skip)
sweep.dz_um     = -3:0.5:3;          % 1-D sweep of defocus   ([] to skip)
sweep.map_dx_um = -2:0.25:2;         % 2-D map: x offsets     ([] to skip the map)
sweep.map_dy_um = -0.8:0.1:0.8;      %          y offsets
sweep.mfd_x_um  = 0.6:0.2:3.0;       % 2-D map over Gaussian mode field diameter (1/e^2 intensity, um)
sweep.mfd_y_um  = 0.3:0.1:1.5;       %   x / y  ([] to skip; always uses a Gaussian spot, whatever cfg.spot_source is)
sweep.base_dz   = 0;                 % defocus held fixed for the offset sweeps / map
sweep.n_lines   = 6;                 % number of strongest modes drawn in the line plots

if exist('cfg_override', 'var'),   cfg = modal_default_cfg(cfg_override); end
if exist('sweep_override', 'var')
    fn = fieldnames(sweep_override);
    for i = 1:numel(fn), sweep.(fn{i}) = sweep_override.(fn{i}); end
end
um = char(181);

%% ---------------- LOAD ----------------
P = load_modal_problem(cfg);
if isempty(P), return; end
nM = P.nM;
base_shift = cfg.beam_offset_um;
if cfg.align_centroids, base_shift = base_shift + P.align_shift; end

% Precompute trapezoid weights and the mode norms once; each sweep point is
% then a single matrix-vector product.
nx = numel(P.intX);  ny = numel(P.intY);
wx = ones(1, nx) * (P.intX(2) - P.intX(1));  wx([1 end]) = wx([1 end]) / 2;
wy = ones(ny, 1) * (P.intY(2) - P.intY(1));  wy([1 end]) = wy([1 end]) / 2;
W  = wy * wx;
M  = reshape(P.E_modes, [], nM);
[~, ~, gram, Pm] = mode_overlap_fractions(P.intX, P.intY, zeros(ny, nx), P.E_modes, 1, P.E_other);
offd = gram - eye(nM);
fprintf('Mode orthogonality: max off-diagonal Gram = %.3f\n', max(abs(offd(:))));

evalpt = @(offs, dz) mode_fracs(P, M, W, Pm, base_shift + offs, dz, cfg.extra_phase_fn);

%% ---------------- BASE POINT ----------------
f0 = evalpt([0 0], sweep.base_dz);
fprintf('Base point (offset 0,0, dz=%.2f): TM00-ish mode 1 = %.2f %%, other listed = %.2f %%, unaccounted = %.2f %%\n', ...
        sweep.base_dz, 100 * f0(1), 100 * sum(f0(2:end)), 100 * max(0, 1 - sum(f0)));

%% ---------------- 1-D SWEEPS ----------------
sw = struct('name', {'x offset', 'y offset', 'defocus'}, ...
            'vals', {sweep.dx_um, sweep.dy_um, sweep.dz_um}, 'F', []);
for s = 1:3
    v = sw(s).vals;  if isempty(v), continue; end
    F = zeros(nM, numel(v));
    for i = 1:numel(v)
        switch s
            case 1, F(:, i) = evalpt([v(i) 0], sweep.base_dz);
            case 2, F(:, i) = evalpt([0 v(i)], sweep.base_dz);
            case 3, F(:, i) = evalpt([0 0], v(i));
        end
    end
    sw(s).F = F;
    fprintf('Sweep %-9s done (%d points)\n', sw(s).name, numel(v));
end

%% ---------------- 2-D OFFSET MAP ----------------
doMap = ~isempty(sweep.map_dx_um) && ~isempty(sweep.map_dy_um);
if doMap
    mx = sweep.map_dx_um;  my = sweep.map_dy_um;
    F1 = zeros(numel(my), numel(mx));  FH = F1;  FU = F1;
    for iy = 1:numel(my)
        for ix = 1:numel(mx)
            f = evalpt([mx(ix) my(iy)], sweep.base_dz);
            F1(iy, ix) = f(1);  FH(iy, ix) = sum(f(2:end));  FU(iy, ix) = max(0, 1 - sum(f));
        end
        fprintf('  map row %d/%d\n', iy, numel(my));
    end
end

%% ---------------- MODE FIELD DIAMETER MAP (Gaussian spot) ----------------
doMfd = ~isempty(sweep.mfd_x_um) && ~isempty(sweep.mfd_y_um);
if doMfd
    gx = sweep.mfd_x_um;  gy = sweep.mfd_y_um;
    G1 = zeros(numel(gy), numel(gx));  GH = G1;  GU = G1;
    for iy = 1:numel(gy)
        for ix = 1:numel(gx)
            cfgG = cfg;  cfgG.gauss.mfd_x_um = gx(ix);  cfgG.gauss.mfd_y_um = gy(iy);
            P2 = modal_set_spot(P, make_gaussian_spot(cfgG));
            sh = cfg.beam_offset_um;  if cfg.align_centroids, sh = sh + P2.align_shift; end
            f = mode_fracs(P2, M, W, Pm, sh, sweep.base_dz, cfg.extra_phase_fn);
            G1(iy, ix) = f(1);  GH(iy, ix) = sum(f(2:end));  GU(iy, ix) = max(0, 1 - sum(f));
        end
    end
    [bestG, ib] = max(G1(:));  [iyb, ixb] = ind2sub(size(G1), ib);
    fprintf('MFD map: best coupling into the fundamental = %.1f %% at MFD x = %.2f, y = %.2f %sm (of %d x %d points)\n', ...
            100 * bestG, gx(ixb), gy(iyb), um, numel(gx), numel(gy));
end

%% ---------------- PLOTS ----------------
tags = P.modeTag;
for k = 1:nM
    tags{k} = sprintf('m%02d %s', P.modeIdx(k), tags{k});
end

% choose the modes to draw: those with the largest peak over all 1-D sweeps
peak = zeros(nM, 1);
for s = 1:3
    if ~isempty(sw(s).F), peak = max(peak, max(sw(s).F, [], 2)); end
end
[~, ord] = sort(peak, 'descend');
show = sort(ord(1:min(sweep.n_lines, nM)));
col = lines(numel(show));

figure('Color', 'w', 'Position', [40 40 1500 460], 'Name', 'Sweeps');
units = {[um 'm'], [um 'm'], [um 'm']};
for s = 1:3
    subplot(1, 3, s);  hold on;  box on;  grid on;
    if isempty(sw(s).F), title([sw(s).name ' (skipped)']); continue; end
    v = sw(s).vals;
    for j = 1:numel(show)
        plot(v, 100 * sw(s).F(show(j), :), '-o', 'Color', col(j, :), 'LineWidth', 1.4, ...
             'MarkerSize', 3, 'DisplayName', tags{show(j)});
    end
    plot(v, 100 * sum(sw(s).F, 1), 'k--', 'LineWidth', 1.4, 'DisplayName', 'sum of listed');
    plot(v, 100 * sum(sw(s).F(2:end, :), 1), 'k:', 'LineWidth', 1.8, 'DisplayName', 'all higher-order');
    xlabel([sw(s).name ' (' units{s} ')']);  ylabel('Fraction of incident power (%)');
    title(sw(s).name);
    if s == 1, legend('Location', 'best', 'FontSize', 7); end
end
sgtitle(sprintf('%s   |   %s', P.fileI, P.modeSetName), 'Interpreter', 'none', 'FontSize', 9);

if doMap
    figure('Color', 'w', 'Position', [60 60 1500 420], 'Name', 'Offset map');
    data = {100 * F1, 100 * FH, 100 * FU};
    ttl  = {'Fundamental mode (m01) (%)', 'All higher-order modes (%)', 'Unaccounted (%)'};
    for s = 1:3
        subplot(1, 3, s);
        imagesc(mx, my, data{s}); set(gca, 'YDir', 'normal'); colorbar; colormap(gca, parula);
        xlabel(['x offset (' um 'm)']); ylabel(['y offset (' um 'm)']); title(ttl{s});
    end
    sgtitle(sprintf('Beam offset map (dz = %.2f %sm)', sweep.base_dz, um));
end

if doMfd
    figure('Color', 'w', 'Position', [80 80 1500 420], 'Name', 'MFD map');
    data = {100 * G1, 100 * GH, 100 * GU};
    ttl  = {'Fundamental mode (m01) (%)', 'All higher-order modes (%)', 'Unaccounted (%)'};
    for q = 1:3
        subplot(1, 3, q);
        imagesc(gx, gy, data{q}); set(gca, 'YDir', 'normal'); colorbar; colormap(gca, parula);  hold on;
        plot(gx(ixb), gy(iyb), 'w*', 'MarkerSize', 9);
        xlabel(['MFD x (' um 'm)']); ylabel(['MFD y (' um 'm)']); title(ttl{q});
    end
    sgtitle(sprintf('Gaussian spot, power vs mode field diameter  (white * = best fundamental coupling, %.1f %%)', 100 * bestG));
end

%% ---------------- SAVE ----------------
if ~isempty(cfg.save_figs_dir)
    figs = findall(0, 'Type', 'figure');
    for fg = figs'
        saveas(fg, fullfile(cfg.save_figs_dir, [regexprep(get(fg, 'Name'), '\W+', '_') '.png']));
    end
end

%% ---------------- local functions ----------------
function f = mode_fracs(P, M, W, Pm, shift, dz, extra_phase_fn)
% Fraction of the incident power in each mode for one beam position / defocus.
    E_in = modal_input_field(P, shift, dz, extra_phase_fn);
    ov = M' * (W(:) .* E_in(:));                 % <E_in, m_k> for all modes
    f = abs(ov).^2 ./ (P.P_in_total * Pm);
end

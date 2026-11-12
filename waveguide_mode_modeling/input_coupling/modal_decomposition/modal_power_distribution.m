%% Modal Power Distribution of an Input Spot
% How is the energy of an input irradiance spot (optionally with a phase front)
% distributed across a directory of waveguide modes for ONE geometry?
%
% Inputs
%   1. Zemax POP irradiance file (.txt)           (required)
%   1b. Zemax POP phase file (.txt)               (optional; none = flat phase)
%   2. COMSOL modes for ONE geometry, either
%      (a) a single combined Export>Data .csv: X, Y, then per mode
%          real(ewfd.Ex), imag(ewfd.Ex), real(ewfd.Ey), imag(ewfd.Ey)
%          (neff is read from the column headers; modes are sorted by
%          descending neff so mode 01 = fundamental), or
%      (b) a directory of per-mode files named mode<NN>_<tag>_neff<x.xxxx>.csv
%          e.g. mode01_tm00_neff2.0123.csv, columns x,y,Re/Im(Ex),Re/Im(Ey)
%
% (Alignment / defocus sweeps: see modal_power_sweep.m.)
% Output: bar chart of fraction of incident power per mode (+ "unaccounted"),
% mode-orthogonality (Gram) check, and mode montage.
% Compatible with MATLAB R2018b.

clearvars -except cfg_override; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'helpers', 'waveguide'));

%% ---------------- USER SETTINGS (defaults in modal_default_cfg.m) ----------------
cfg = modal_default_cfg();
% e.g.  cfg.beam_offset_um = [0.5 0];   cfg.defocus_um = 3;   cfg.pol_component = 'Ey';
if exist('cfg_override', 'var')   % optional: struct of settings set before running (batch use)
    cfg = modal_default_cfg(cfg_override);
end
um = char(181);                   % micro sign, avoids file-encoding problems

%% ---------------- 1-3. LOAD INPUT SPOT, MODES, COMMON GRID ----------------
P = load_modal_problem(cfg);
if isempty(P), return; end
intX = P.intX;  intY = P.intY;  nM = P.nM;
E_modes = P.E_modes;  E_other = P.E_other;
neff = P.neff;  modeIdx = P.modeIdx;  modeTag = P.modeTag;  polFrac = P.polFrac;  pol = P.pol;
P_in_total = P.P_in_total;  fileI = P.fileI;

%% ---------------- 4. PLACE INPUT FIELD ON GRID ----------------
shift = cfg.beam_offset_um;
if cfg.align_centroids, shift = shift + P.align_shift; end
fprintf('Beam shift applied: (%.2f, %.2f) %sm,  defocus %.2f %sm\n', shift(1), shift(2), um, cfg.defocus_um, um);
E_in = modal_input_field(P, shift, cfg.defocus_um, cfg.extra_phase_fn);

%% ---------------- 5. OVERLAPS ----------------
[frac, ~, gram, ~] = mode_overlap_fractions(intX, intY, E_in, E_modes, P_in_total, E_other);
P_window = trapz(intY, trapz(intX, abs(E_in).^2, 2));
clipFrac = 1 - P_window / P_in_total;
unaccounted = max(0, 1 - sum(frac));

fprintf('\n%-6s %-6s %-8s %-9s %-10s\n', 'Mode', 'Tag', 'neff', 'pol(%)', 'Power(%)');
for k = 1:nM
    fprintf('%-6d %-6s %-8.4f %-9.1f %-10.3f\n', modeIdx(k), modeTag{k}, neff(k), 100 * polFrac(k), 100 * frac(k));
end
fprintf('Sum of listed modes : %.2f %%\n', 100 * sum(frac));
fprintf('Unaccounted         : %.2f %%  (input outside COMSOL window: %.2f %%)\n', ...
        100 * unaccounted, 100 * clipFrac);
offdiag = gram - eye(nM);
fprintf('Max mode non-orthogonality (off-diagonal Gram): %.3f\n', max(abs(offdiag(:))));
if max(abs(offdiag(:))) > 0.05
    fprintf(2, 'WARNING: modes are not orthogonal -- fractions are not strictly additive.\n');
end
if sum(frac) > 1.02
    fprintf(2, 'WARNING: fractions sum > 100%% (non-orthogonal set or grid problem).\n');
end

%% ---------------- 6. PLOTS ----------------
% Fig 1: input field
figure('Color', 'w', 'Position', [50 50 1100 420], 'Name', 'Input spot');
subplot(1, 2, 1);
imagesc(intX, intY, abs(E_in).^2); axis image; set(gca, 'YDir', 'normal'); colorbar; colormap(gca, jet);
xlabel(['X (' um 'm)']); ylabel(['Y (' um 'm)']); title('Input irradiance on window'); hold on;
contour(intX, intY, abs(E_modes(:, :, 1)).^2, 3, 'w');
subplot(1, 2, 2);
ph = angle(E_in);  ph(abs(E_in).^2 < 0.01 * max(abs(E_in(:)).^2)) = NaN;
imagesc(intX, intY, ph, 'AlphaData', ~isnan(ph)); axis image; set(gca, 'YDir', 'normal');
colorbar; colormap(gca, hsv); caxis([-pi pi]);
xlabel(['X (' um 'm)']); ylabel(['Y (' um 'm)']); title('Input phase (rad), >1% intensity');

% Fig 2: bar chart
labels = cell(nM + 1, 1);
for k = 1:nM
    if isempty(modeTag{k}), nameStr = sprintf('m%02d', modeIdx(k)); else, nameStr = modeTag{k}; end
    if isnan(neff(k)), labels{k} = nameStr;
    else, labels{k} = sprintf('%s (%.4f)', nameStr, neff(k)); end
end
labels{nM + 1} = sprintf('unacc.');
vals = [frac; unaccounted] * 100;

figure('Color', 'w', 'Position', [100 100 1100 600], 'Name', 'Modal power distribution');
for p = 1:2
    subplot(2, 1, p);
    b = bar(1:nM + 1, vals, 0.7, 'FaceColor', 'flat', 'EdgeColor', 'none');
    b.CData = repmat([0.15 0.45 0.75], nM + 1, 1);
    b.CData(end, :) = [0.65 0.65 0.65];
    set(gca, 'XTick', 1:nM + 1, 'XTickLabel', labels, 'XTickLabelRotation', 60, 'FontSize', 8);
    ylabel('Fraction of incident power (%)'); grid on; box on;
    if p == 1
        title(sprintf('%s  |  sum of modes = %.1f%%, unaccounted = %.1f%%  |  pol = %s', ...
              fileI, 100 * sum(frac), 100 * unaccounted, pol), 'Interpreter', 'none');
        for k = 1:nM + 1
            text(k, vals(k), sprintf('%.2f', vals(k)), 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', 'FontSize', 8);
        end
    else
        set(gca, 'YScale', 'log'); ylim([1e-4 100]);
        title('Same data, log scale (small higher-order contributions)');
    end
end

% Fig 3: orthogonality + mode montage
nShow = min(nM, cfg.max_montage_modes);
nc = ceil(sqrt(nShow + 1));  nr = ceil((nShow + 1) / nc);
figure('Color', 'w', 'Position', [150 150 1100 750], 'Name', 'Modes and orthogonality');
subplot(nr, nc, 1);
imagesc(gram); axis image; colorbar; caxis([0 1]); colormap(gca, parula);
title('|<mi,mj>| (Gram)'); xlabel('mode #'); ylabel('mode #');
amax = max(abs(E_modes), [], 3) + max(abs(E_other), [], 3);
amax = amax / max(amax(:));
cols_ = find(any(amax > 0.02, 1));  rows_ = find(any(amax > 0.02, 2));
mx = 0.1 * (intX(cols_(end)) - intX(cols_(1)));  my = 0.3 * (intY(rows_(end)) - intY(rows_(1)));
cropX = [intX(cols_(1)) - mx, intX(cols_(end)) + mx];
cropY = [intY(rows_(1)) - my, intY(rows_(end)) + my];
cmap = [linspace(0, 1, 64)', linspace(0, 1, 64)', ones(64, 1); ...
        ones(64, 1), linspace(1, 0, 64)', linspace(1, 0, 64)'];
for k = 1:nShow
    subplot(nr, nc, k + 1);
    Ek = E_modes(:, :, k);
    [~, imax] = max(abs(Ek(:)));
    Ek = Ek * exp(-1i * angle(Ek(imax)));            % fix arbitrary global phase for display
    imagesc(intX, intY, real(Ek) / max(abs(Ek(:)))); axis image; set(gca, 'YDir', 'normal');
    caxis([-1 1]); colormap(gca, cmap); xlim(cropX); ylim(cropY);
    title(sprintf('%s  %.1f%%', strtrim([sprintf('m%02d ', modeIdx(k)) modeTag{k}]), 100 * frac(k)), 'FontSize', 9);
    set(gca, 'XTick', [], 'YTick', []);
end
sgtitle('Mode fields (Re, signed) and orthogonality check');

%% ---------------- 7. OPTIONAL: SAVE FIGURES ----------------
if ~isempty(cfg.save_figs_dir)
    figs = findall(0, 'Type', 'figure');
    for fg = figs'
        saveas(fg, fullfile(cfg.save_figs_dir, [regexprep(get(fg, 'Name'), '\W+', '_') '.png']));
    end
end

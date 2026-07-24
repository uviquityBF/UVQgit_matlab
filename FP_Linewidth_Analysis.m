function FP_Linewidth_Analysis()
% Laser linewidth measurement from a scanning Fabry-Perot interferometer.
% Reads a single-channel waveform from a Rohde & Schwarz MXO44 CSV export.
%
% Method: FWHM of one FP transmission peak / spacing between two adjacent peaks * FSR
% Instrument: Thorlabs SA200-3B  (FSR = 1.5 GHz)
%
% USAGE:
%   1. Set csv_file to your file path.
%   2. Run. The overview plot shows all detected peaks numbered in time order.
%   3. Set peak_indices to the two adjacent peaks you want to analyze.
%   4. Run again for the detail plot and linewidth result.

    close all; clear; clc;

    %% --- USER SETTINGS ---
    csv_file         = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_07_20_LabData\linewidth measurment\sa200-3B\DLPro_at_250mA\Dlpro-sa2003b-test2_2026-07-21_0_192721.csv'; 
    % full path to the MXO44 CSV file
    FSR_GHz          = 1.5;            % SA200-3B free spectral range (GHz)

    % Peak detection — adjust if peaks are missed or noise is picked up
    min_prominence_V = 0.02;           % min peak height above local background (V)
    min_sep_ms       = 0.5;            % min time between adjacent peaks (ms)

    % Which two peaks to analyze — set after viewing the overview plot
    peak_indices    = [1, 2];          % [left_peak, right_peak] by time order

    % Optional time window to restrict the view (ms). Leave [] for full waveform.
    t_window_ms     = [];              % e.g. [2.0, 8.0]

    %% --- 1. LOAD DATA ---
    fprintf('Loading %s ...\n', csv_file);

    % Pass 1: count lines to the TIME header and detect delimiter
    fid = fopen(csv_file, 'r');
    if fid < 0, error('Cannot open file: %s', csv_file); end
    n_skip   = 0;
    line_str = '';
    while ~feof(fid) && ~is_time_header(line_str)
        line_str = fgetl(fid);
        n_skip   = n_skip + 1;
    end
    first_data_line = fgetl(fid);   % peek at first data line for diagnosis
    fclose(fid);

    if n_skip == 0
        error('Could not find a line starting with TIME in the CSV file.');
    end
    fprintf('  TIME header on line %d\n', n_skip);
    fprintf('  First data line: [%s]\n', first_data_line(1:min(80, length(first_data_line))));

    if contains(line_str, ','),     delim = ',';
    elseif contains(line_str, ';'), delim = ';';
    else,                           delim = '\t';
    end
    fprintf('  Delimiter: "%s"\n', delim);

    % Pass 2: reopen, skip n_skip lines (header + TIME row), read data fresh
    fid = fopen(csv_file, 'r');
    for k = 1:n_skip
        fgetl(fid);
    end
    raw = textscan(fid, '%f %f', 'Delimiter', delim, 'CollectOutput', true);
    fclose(fid);

    fprintf('  Rows parsed: %d\n', size(raw{1}, 1));
    if size(raw{1}, 1) < 2
        error('textscan read fewer than 2 rows. Check the first data line printed above.');
    end

    t_s  = raw{1}(:, 1);
    v    = raw{1}(:, 2);
    dt_s = t_s(2) - t_s(1);
    t_ms = t_s * 1e3;

    fprintf('  Loaded %.2f ms, %.0f MSa/s, %d samples\n', ...
            (t_s(end) - t_s(1)) * 1e3, 1/dt_s/1e6, length(t_s));

    %% --- 2. OPTIONAL TIME WINDOW ---
    if ~isempty(t_window_ms)
        mask = t_ms >= t_window_ms(1) & t_ms <= t_window_ms(2);
        t_ms = t_ms(mask);
        v    = v(mask);
    end

    %% --- 3. DETECT PEAKS ---
    min_sep_samples = round(min_sep_ms * 1e-3 / dt_s);
    [pks, locs] = findpeaks(v, ...
                            'MinPeakProminence', min_prominence_V, ...
                            'MinPeakDistance',   min_sep_samples);

    fprintf('  Found %d peaks (prominence > %.3f V, separation > %.2f ms)\n', ...
            length(pks), min_prominence_V, min_sep_ms);
    if length(pks) < 2
        error('Fewer than 2 peaks found. Lower min_prominence_V or min_sep_ms.');
    end

    %% --- 4. OVERVIEW PLOT (numbered peaks) ---
    % Subsample for fast plotting if the waveform is very long
    stride = max(1, floor(length(t_ms) / 20000));
    figure('Color','w','Name','FP Overview — check peak numbers', ...
           'Position',[50, 450, 1100, 350]);
    plot(t_ms(1:stride:end), v(1:stride:end), 'b', 'LineWidth', 0.6);
    hold on;
    plot(t_ms(locs), pks, 'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
    for p = 1:length(pks)
        text(t_ms(locs(p)), pks(p), sprintf('  %d', p), ...
             'FontSize', 8, 'Color', 'r', 'VerticalAlignment', 'bottom');
    end
    xlabel('Time (ms)'); ylabel('Voltage (V)');
    title(sprintf('FP Waveform Overview — %d peaks detected. Set peak\\_indices then re-run.', length(pks)));
    grid on;

    %% --- 5. SELECT THE TWO ANALYSIS PEAKS ---
    if any(peak_indices > length(pks))
        error('peak_indices [%d %d] out of range — only %d peaks found.', ...
              peak_indices(1), peak_indices(2), length(pks));
    end
    p1 = peak_indices(1);
    p2 = peak_indices(2);

    t_p1 = t_ms(locs(p1));    v_p1 = pks(p1);
    t_p2 = t_ms(locs(p2));    v_p2 = pks(p2);
    spacing_ms = t_p2 - t_p1;

    fprintf('\nAnalysis peaks:\n');
    fprintf('  Peak %d: t = %.5f ms,  V = %.4f V\n', p1, t_p1, v_p1);
    fprintf('  Peak %d: t = %.5f ms,  V = %.4f V\n', p2, t_p2, v_p2);
    fprintf('  Spacing: %.5f ms  (= %.1f GHz)\n', spacing_ms, FSR_GHz);

    %% --- 6. FWHM OF PEAK 1 ---
    half_win_ms = 0.35 * spacing_ms;
    win = t_ms >= (t_p1 - half_win_ms) & t_ms <= (t_p1 + half_win_ms);
    t_win = t_ms(win);
    v_win = v(win);

    [fwhm_ms, t_left_ms, t_right_ms, half_max_V] = compute_fwhm(t_win, v_win);

    if isnan(fwhm_ms)
        warning('FWHM could not be determined for peak %d. Try adjusting the window or peak.', p1);
        linewidth_MHz = NaN;
    else
        linewidth_MHz = (fwhm_ms / spacing_ms) * FSR_GHz * 1e3;
        fprintf('\n=== LINEWIDTH RESULT ===\n');
        fprintf('  FWHM (time):      %.5f ms\n', fwhm_ms);
        fprintf('  Peak spacing:     %.5f ms  (= %.1f GHz FSR)\n', spacing_ms, FSR_GHz);
        fprintf('  Laser linewidth:  %.2f MHz\n', linewidth_MHz);
    end

    %% --- 7. DETAIL PLOT ---
    margin_ms   = 0.6 * spacing_ms;
    det_mask    = t_ms >= (t_p1 - margin_ms) & t_ms <= (t_p2 + margin_ms);
    t_det       = t_ms(det_mask);
    v_det       = v(det_mask);

    % Calibrated frequency axis: peak 1 is origin, peak 2 is +FSR
    GHz_per_ms  = FSR_GHz / spacing_ms;
    freq_GHz    = (t_det - t_p1) * GHz_per_ms;

    figure('Color','w','Name','FP Peak Analysis','Position',[100, 50, 950, 520]);
    plot(freq_GHz, v_det, 'b', 'LineWidth', 1.2);
    hold on;

    % Mark the two selected peaks
    plot(0,                      v_p1, 'rv', 'MarkerFaceColor','r','MarkerSize',11);
    plot(FSR_GHz,                v_p2, 'rv', 'MarkerFaceColor','r','MarkerSize',11);

    % FSR reference arrow
    y_fsr  = min(v_det) + 0.08 * (max(v_det) - min(v_det));
    line([0, FSR_GHz], [y_fsr, y_fsr], 'Color','k','LineWidth',1.5,'LineStyle','-');
    plot(0,       y_fsr, 'k>', 'MarkerFaceColor','k','MarkerSize',6);
    plot(FSR_GHz, y_fsr, 'k<', 'MarkerFaceColor','k','MarkerSize',6);
    text(FSR_GHz/2, y_fsr, sprintf('  FSR = %.1f GHz  ', FSR_GHz), ...
         'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'FontWeight','bold');

    % FWHM bracket on peak 1
    if ~isnan(fwhm_ms)
        fwhm_GHz     = fwhm_ms * GHz_per_ms;
        t_left_GHz   = (t_left_ms  - t_p1) * GHz_per_ms;
        t_right_GHz  = (t_right_ms - t_p1) * GHz_per_ms;
        line([t_left_GHz, t_right_GHz], [half_max_V, half_max_V], ...
             'Color','m','LineWidth',2.5);
        tick_h = 0.03 * (max(v_det) - min(v_det));
        line([t_left_GHz,  t_left_GHz],  [half_max_V - tick_h, half_max_V + tick_h], 'Color','m','LineWidth',2);
        line([t_right_GHz, t_right_GHz], [half_max_V - tick_h, half_max_V + tick_h], 'Color','m','LineWidth',2);
        text(t_left_GHz, half_max_V, sprintf('FWHM = %.1f MHz  ', linewidth_MHz), ...
             'Color','m','FontSize',11,'FontWeight','bold', ...
             'HorizontalAlignment','right','VerticalAlignment','middle');
    end

    xlabel('Relative Frequency (GHz)','FontSize',12);
    ylabel('Voltage (V)','FontSize',12);
    if ~isnan(fwhm_ms)
        title(sprintf('FP Transmission — Laser Linewidth = %.2f MHz  (peaks %d & %d)', ...
                      linewidth_MHz, p1, p2), 'FontSize',12);
    else
        title(sprintf('FP Transmission — FWHM not determined  (peaks %d & %d)', p1, p2));
    end
    grid on;
end


%% =========================================================================
%  Returns true only if line_str is the TIME column header row
%  (i.e. starts with 'TIME' followed immediately by a delimiter, not more letters).
%  This prevents matching 'TimebaseScale' or similar metadata keys.
% =========================================================================
function result = is_time_header(line_str)
    if ~ischar(line_str) || length(line_str) < 4
        result = false;
        return;
    end
    s = strtrim(line_str);
    result = strcmpi(s(1:4), 'TIME') && (length(s) == 4 || ~isletter(s(5)));
end

%% =========================================================================
%  FWHM — finds half-maximum crossings by linear interpolation.
%  Baseline is estimated as the minimum of the window (works for sharp peaks
%  on a smooth low background). Returns NaN if crossings cannot be found.
% =========================================================================
function [fwhm, t_left, t_right, half_max] = compute_fwhm(t, v)
    fwhm = NaN; t_left = NaN; t_right = NaN; half_max = NaN;

    [v_peak, idx_peak] = max(v);
    baseline  = min(v);
    half_max  = baseline + (v_peak - baseline) / 2;

    % Left crossing (search left from peak)
    t_left = NaN;
    for i = idx_peak:-1:2
        if v(i-1) <= half_max && v(i) >= half_max
            t_left = interp1([v(i-1), v(i)], [t(i-1), t(i)], half_max);
            break;
        end
    end

    % Right crossing (search right from peak)
    t_right = NaN;
    for i = idx_peak:length(v)-1
        if v(i) >= half_max && v(i+1) <= half_max
            t_right = interp1([v(i), v(i+1)], [t(i), t(i+1)], half_max);
            break;
        end
    end

    if ~isnan(t_left) && ~isnan(t_right)
        fwhm = t_right - t_left;
    end
end

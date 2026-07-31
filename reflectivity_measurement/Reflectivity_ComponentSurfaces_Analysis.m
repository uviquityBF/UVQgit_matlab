function Reflectivity_ComponentSurfaces_Analysis()
% Single-bounce surface reflectivity from Ocean Optics (Maya) spectra.
%
% For each source (deuterium lamp / blue laser), takes the ratio of a
% single-bounce-off-surface spectrum ("posN") to a no-bounce reference
% spectrum ("RefThru") captured through air. Each reflected spectrum is
% paired with whichever reference was closest in time, since the recorded
% lamp shape drifts a bit with collection angle/alignment.
%
% Because that drift means the raw ratio isn't a clean flat reflectivity
% curve, each ratio is renormalized to its own mean value over a visible
% wavelength band (vis_norm_band_nm). This gives a *relative* reflectivity:
% how much lower the UV reflects compared to the visible, for that surface.
% If you have an absolute reflectivity number for the visible band from
% another measurement, set R_vis_absolute_pct to convert to an absolute
% UV reflectivity estimate.
%
% USAGE:
%   1. Set data_dir to the folder of spectra.
%   2. Run. Check the raw-spectra figure to confirm sensible source/role
%      grouping, and the summary table printed to the console.
%   3. Adjust vis_norm_band_nm / uv_report_band_nm as needed and re-run.

    close all; clear; clc;

    %% --- USER SETTINGS ---
    data_dir = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_07_27_LabData\reflectivity_component_surfaces';

    vis_norm_band_nm    = [450, 650];  % band used to normalize each ratio curve to 1 (relative reflectivity anchor)
    uv_report_band_nm   = [200, 250];  % band of interest for the UV reflectivity summary
    plot_range_nm       = [180, 900];  % x-axis range for the raw-spectra and ratio plots
    smooth_window_nm    = 2;           % moving-average smoothing applied before taking ratios
    min_valid_signal_frac = 0.02;      % mask the ratio where the reference signal is below this fraction of its own peak

    R_vis_absolute_pct  = NaN;         % known absolute reflectivity (%) in vis_norm_band_nm from another measurement.
                                        % Leave NaN to report relative-only. The blue laser peak ratio (printed in the
                                        % summary) is a candidate for this number if its wavelength falls in vis_norm_band_nm.

    %% --- 1. LOAD ALL SPECTRA ---
    files = dir(fullfile(data_dir, '*.txt'));
    if isempty(files)
        error('No .txt files found in %s', data_dir);
    end
    fprintf('Loading %d spectra from %s ...\n', numel(files), data_dir);

    S = repmat(struct('name','', 'wl',[], 'counts',[], 'Tint',NaN, 'Nscans',NaN, ...
                       'source','', 'role','', 'pos','', 't_s',NaN, 'cps',[], 'cps_smooth',[]), 1, numel(files));
    for k = 1:numel(files)
        DAT = readOceanOptics_SpectrumFile(fullfile(data_dir, files(k).name));
        meta = parse_filename(files(k).name);

        S(k).name   = files(k).name;
        S(k).wl     = DAT.wl;
        S(k).counts = DAT.counts;
        S(k).Tint   = DAT.Tint;
        S(k).Nscans = DAT.Nscans;
        S(k).source = meta.source;
        S(k).role   = meta.role;
        S(k).pos    = meta.pos;
        S(k).t_s    = meta.t_s;

        fprintf('  %-70s  source=%-10s role=%-10s pos=%-6s\n', files(k).name, meta.source, meta.role, meta.pos);
    end

    %% --- 2. COMMON WAVELENGTH GRID ---
    wl_ref = S(1).wl;
    for k = 2:numel(S)
        if numel(S(k).wl) ~= numel(wl_ref) || max(abs(S(k).wl - wl_ref)) > 1e-3
            warning('%s wavelength grid differs from the first file -- interpolating onto a common grid.', S(k).name);
            S(k).counts = interp1(S(k).wl, S(k).counts, wl_ref, 'linear', 'extrap');
            S(k).wl = wl_ref;
        end
    end

    %% --- 3. NORMALIZE TO COUNTS/SEC AND SMOOTH ---
    dwl = median(diff(wl_ref));
    win_pts = max(1, round(smooth_window_nm / dwl));
    for k = 1:numel(S)
        S(k).cps = S(k).counts / (S(k).Tint * S(k).Nscans);
        S(k).cps_smooth = movmean(S(k).cps, win_pts);
    end

    %% --- 4. RATIO EACH REFLECTED SPECTRUM AGAINST ITS NEAREST-TIME REFERENCE ---
    sources = unique({S.source});
    results = struct('source',{}, 'pos',{}, 'filename',{}, 'ref_filename',{}, ...
                      'wl',{}, 'ratio',{}, 'ratio_norm',{}, 'peak_wl',{}, 'peak_ratio',{}, ...
                      'uv_mean_relative',{}, 'uv_mean_absolute_pct',{});

    for si = 1:numel(sources)
        src = sources{si};
        idx_ref  = find(strcmp({S.source}, src) & strcmp({S.role}, 'reference'));
        idx_refl = find(strcmp({S.source}, src) & strcmp({S.role}, 'reflected'));
        if isempty(idx_ref) || isempty(idx_refl)
            warning('Source "%s" is missing a reference or a reflected spectrum -- skipping.', src);
            continue;
        end

        for jj = 1:numel(idx_refl)
            k = idx_refl(jj);
            [~, bi] = min(abs([S(idx_ref).t_s] - S(k).t_s));
            kref = idx_ref(bi);

            ref_sig  = S(kref).cps_smooth;
            refl_sig = S(k).cps_smooth;
            valid = abs(ref_sig) > min_valid_signal_frac * max(abs(ref_sig));

            ratio = nan(size(wl_ref));
            ratio(valid) = refl_sig(valid) ./ ref_sig(valid);

            vis_mask = wl_ref >= vis_norm_band_nm(1) & wl_ref <= vis_norm_band_nm(2) & valid;
            norm_factor = mean(ratio(vis_mask), 'omitnan');
            ratio_norm = ratio / norm_factor;

            [~, pk] = max(refl_sig);
            peak_wl = wl_ref(pk);
            peak_mask = abs(wl_ref - peak_wl) <= 3 & valid;
            peak_ratio = mean(refl_sig(peak_mask)) / mean(ref_sig(peak_mask));

            uv_mask = wl_ref >= uv_report_band_nm(1) & wl_ref <= uv_report_band_nm(2) & valid;
            uv_mean_relative = mean(ratio_norm(uv_mask), 'omitnan');
            if isnan(R_vis_absolute_pct)
                uv_mean_absolute_pct = NaN;
            else
                uv_mean_absolute_pct = uv_mean_relative * R_vis_absolute_pct;
            end

            results(end+1) = struct('source', src, 'pos', S(k).pos, 'filename', S(k).name, ...
                'ref_filename', S(kref).name, 'wl', wl_ref, 'ratio', ratio, 'ratio_norm', ratio_norm, ...
                'peak_wl', peak_wl, 'peak_ratio', peak_ratio, ...
                'uv_mean_relative', uv_mean_relative, 'uv_mean_absolute_pct', uv_mean_absolute_pct); %#ok<AGROW>
        end
    end

    %% --- 5. SUMMARY TABLE ---
    fprintf('\n=== Reflectivity summary ===\n');
    fprintf('Relative UV/vis is the UV-band ratio after normalizing to the %g-%g nm band (i.e. UV reflectivity relative to visible).\n', vis_norm_band_nm);
    fprintf('Peak ratio is the raw (unnormalized) reflected/reference ratio at each spectrum''s own peak wavelength -- meaningful mainly for the narrowband laser.\n\n');
    fprintf('%-11s %-8s %6s %10s %8s %10s   %s\n', 'Source', 'Pos', 'Pk[nm]', 'PeakRatio', 'UV/Vis', 'AbsUV[%]', 'Reference used');
    for ir = 1:numel(results)
        r = results(ir);
        if isnan(r.uv_mean_relative)
            uvstr = 'n/a';   % narrowband source has no signal in the vis-norm band
        else
            uvstr = sprintf('%.3f', r.uv_mean_relative);
        end
        if isnan(r.uv_mean_absolute_pct)
            absstr = '--';
        else
            absstr = sprintf('%.2f', r.uv_mean_absolute_pct);
        end
        fprintf('%-11s %-8s %6.1f %10.3f %8s %10s   %s\n', ...
            r.source, r.pos, r.peak_wl, r.peak_ratio, uvstr, absstr, r.ref_filename);
    end
    fprintf('\n');

    %% --- 6. PLOTS ---
    plot_raw_spectra(S, sources, plot_range_nm);
    plot_ratio(results, sources, plot_range_nm, vis_norm_band_nm, 'ratio', 'Reflected / Reference (raw)', 'Raw Ratio');
    plot_ratio(results, sources, [plot_range_nm(1), vis_norm_band_nm(2)], uv_report_band_nm, 'ratio_norm', ...
        'Reflectivity relative to visible band (=1)', 'Normalized Ratio (UV region)');

end


%% ================= Helper functions =================

function DAT = readOceanOptics_SpectrumFile(filepath)
% Reads an Ocean Optics (SpectraSuite/OceanView) plain-text spectrum file.
    fid = fopen(filepath, 'r');
    if fid < 0
        error('Cannot open file: %s', filepath);
    end

    DAT.Tint = NaN;
    DAT.Nscans = NaN;
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        if contains(line, 'Integration Time')
            tok = regexp(line, '([\d.]+[Ee]?[+-]?\d*)\s*$', 'match', 'once');
            if ~isempty(tok), DAT.Tint = str2double(tok); end
        elseif contains(line, 'Scans to average')
            tok = regexp(line, '(\d+)\s*$', 'match', 'once');
            if ~isempty(tok), DAT.Nscans = str2double(tok); end
        elseif contains(line, 'Begin Spectral Data')
            break;
        end
    end

    raw = textscan(fid, '%f %f');
    fclose(fid);

    DAT.wl = raw{1};
    DAT.counts = raw{2};
end


function meta = parse_filename(fname)
% Extracts source, role (reference/reflected), position label, and capture
% time-of-day (seconds) from the lab's spectrum filename convention.
    if contains(fname, 'DeutLamp', 'IgnoreCase', true)
        meta.source = 'Deuterium';
    elseif contains(fname, 'Nichia', 'IgnoreCase', true)
        meta.source = 'BlueLaser';
    else
        meta.source = 'Unknown';
    end

    if contains(fname, 'RefThru', 'IgnoreCase', true)
        meta.role = 'reference';
    else
        meta.role = 'reflected';
    end

    postok = regexpi(fname, 'pos\d+[a-z]*', 'match', 'once');
    if isempty(postok)
        postok = 'unk';
    end
    meta.pos = postok;

    ttok = regexp(fname, '__\d+__(\d{2})-(\d{2})-(\d{2})-\d+', 'tokens', 'once');
    if ~isempty(ttok)
        hms = str2double(ttok);
        meta.t_s = hms(1)*3600 + hms(2)*60 + hms(3);
    else
        meta.t_s = NaN;
    end
end


function plot_raw_spectra(S, sources, plot_range_nm)
    figure('Name', 'Raw Spectra');
    for si = 1:numel(sources)
        subplot(numel(sources), 1, si); hold on;
        idx = find(strcmp({S.source}, sources{si}));
        for k = idx
            if strcmp(S(k).role, 'reference')
                ls = '--'; lw = 2;
            else
                ls = '-'; lw = 1.2;
            end
            plot(S(k).wl, S(k).cps_smooth, ls, 'LineWidth', lw, ...
                'DisplayName', sprintf('%s (%s)', S(k).pos, S(k).role));
        end
        set(gca, 'YScale', 'log');
        xlim(plot_range_nm); grid on;
        xlabel('Wavelength [nm]'); ylabel('Counts / sec');
        title(sources{si});
        legend('show', 'Location', 'best');
    end
end


function plot_ratio(results, sources, plot_range_nm, highlight_band_nm, field, ylab, fig_name)
    figure('Name', fig_name);
    for si = 1:numel(sources)
        subplot(numel(sources), 1, si); hold on;
        idx = find(strcmp({results.source}, sources{si}));
        if isempty(idx)
            title(sprintf('%s (no data)', sources{si}));
            continue;
        end
        yl = ylim; %#ok<NASGU>
        patch([highlight_band_nm(1), highlight_band_nm(2), highlight_band_nm(2), highlight_band_nm(1)], ...
              [-10, -10, 10, 10], [1 1 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
        yline(1, 'k:', 'HandleVisibility', 'off');
        for ii = idx
            plot(results(ii).wl, results(ii).(field), 'LineWidth', 1.3, 'DisplayName', results(ii).pos);
        end
        xlim(plot_range_nm); ylim([-0.2, 1.5]); grid on;
        xlabel('Wavelength [nm]'); ylabel(ylab);
        title(sources{si});
        legend('show', 'Location', 'best');
    end
end

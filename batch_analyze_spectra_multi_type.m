% %% Open Avantes Spectra (ASCII) with Gaussian and Centroid Analysis  with Smoothing
%% Open Avantes/OceanOptics/CSV Spectra with Gaussian, Centroid, and 3D Plotting
function batch_analyze_spectra_multi_type()
    close all
    % --- Setup Paths ---
    default_dir = 'C:\Users\brent\Downloads\NichiaGC1_20C_450mA_4hr_hold05212026_VBG_spectra';
    if ~exist(default_dir, 'dir'); default_dir = pwd; end
    
    % --- CONFIGURATION ---
    % Toggle between: 'Avantes', 'OceanOptics', or 'CSV_2Col'
    type = 'CSV_2Col'; 
    
    % Update file filter extension based on selection
    if strcmp(type, 'CSV_2Col')
        file_filter = {'*.csv;*.txt', 'CSV/Text Data Files (*.csv, *.txt)'};
    else
        file_filter = {'*.txt', 'Spectra Files (*.txt)'};
    end
    
    % 1. Multi-select files
    cd(default_dir);
    [filenames, start_path] = uigetfile(file_filter, ...
        'Select Spectra Files', default_dir, 'MultiSelect', 'on');
    
    if isequal(filenames, 0); return; end
    if ischar(filenames); filenames = {filenames}; end 
    
    smooth_ans = questdlg('Apply low-pass smoothing (Savitzky-Golay)?', ...
                          'Filter Settings', 'Yes', 'No', 'No');
    do_smooth = strcmp(smooth_ans, 'Yes');

    num_files = length(filenames);
    target_int_time = 100; % Target in ms
    
    % Initialize containers
    all_WL = []; all_normalized = []; all_time_scaled = [];
    labels = cell(1, num_files);
    fit_center_wl = zeros(1, num_files);
    centroid_wl = zeros(1, num_files); 

    fprintf('\n--- Processing %s Data ---\n', type);
    
    for k = 1:num_files
        full_path = fullfile(start_path, filenames{k});
        labels{k} = strrep(filenames{k}, '_', ' ');

        % --- Adaptive Integration Time Extraction ---
        int_time = 100; % Default to target for raw CSV to keep scale factor 1
        
        if ~strcmp(type, 'CSV_2Col')
            fid = fopen(full_path, 'r');
            search_limit = ifempty(strcmp(type, 'Avantes'), 9, 17);
            for i = 1:search_limit
                line_str = fgetl(fid);
                if contains(lower(line_str), 'integration time')
                    if strcmp(type, 'Avantes')
                        nums = sscanf(line_str, '%*[^0-9.]%f');
                        if ~isempty(nums); int_time = nums(1); end 
                    else % OceanOptics
                        nums = sscanf(line_str, 'Integration Time (sec): %e');
                        if ~isempty(nums); int_time = nums(1) * 1000; end 
                    end
                end
            end
            fclose(fid);
        end

        % --- Read Numeric Data ---
        switch type
            % Replaced legacy readmatrix check with dlmread/try-catch for pure backward safety
            case 'Avantes'
                data = dlmread(full_path, ';', 9, 0);
            case 'OceanOptics'
                raw = importdata(full_path, '\t', 17);
                if isstruct(raw), data = raw.data; else, data = raw; end
            case 'CSV_2Col'
                % Directly read 2-column comma-separated or tab-separated file without header parsing
                try
                    data = readmatrix(full_path);
                catch
                    data = dlmread(full_path, ',', 1, 0); % Fallback for older MATLAB versions
                end
        end
        
        WL = data(:,1);
        Intensity = data(:,2);
        
        if do_smooth
            Intensity = sgolayfilt(Intensity, 3, 11);
        end

        if k == 1; all_WL = WL; end 

        % --- Gaussian Fit (430-470nm) ---
        idx = find(WL >= 430 & WL <= 470);
        if ~isempty(idx) && length(idx) > 3
            x = WL(idx); y = Intensity(idx);
            y_fit = y; y_fit(y_fit <= 0) = min(y_fit(y_fit > 0)); 
            logY = log(y_fit);
            coeffs = polyfit(x, logY, 2);
            fit_center_wl(k) = -coeffs(2) / (2 * coeffs(1));
        else
            fit_center_wl(k) = NaN;
        end

        % --- Centroid Calculation ---
        baseMask = (WL >= 200) & (WL <= 300);
        avg_noise = ifempty(any(baseMask), mean(Intensity(baseMask)), 0);
        auto_thresh = avg_noise + (10 * ifempty(any(baseMask), std(Intensity(baseMask)), 10));

        if ~isempty(idx)
            wl_seg = WL(idx); int_seg = Intensity(idx);
            [maxVal, ~] = max(int_seg);
            if maxVal > auto_thresh
                cutOff = avg_noise + (maxVal - avg_noise) * 0.5;
                pMask = int_seg >= cutOff;
                centroid_wl(k) = sum(wl_seg(pMask) .* int_seg(pMask)) / sum(int_seg(pMask));
            else
                centroid_wl(k) = NaN;
            end
        else
            centroid_wl(k) = NaN;
        end
        
        fprintf('%-40s | G: %-7.2f | C: %-7.2f\n', filenames{k}, fit_center_wl(k), centroid_wl(k));

        % Rescaling
        max_val = max(Intensity); if max_val <= 0; max_val = 1; end 
        all_normalized(:, k) = Intensity / max_val;
        all_time_scaled(:, k) = Intensity * (target_int_time / int_time);
    end

    % --- Visualization & Export ---
    plot_spectra(all_WL, all_normalized, ['Normalized (', type, ')'], labels);
    plot_spectra(all_WL, all_time_scaled, sprintf('Time-Corrected (%gms)', target_int_time), labels);
    
    % New 3D Plot Layout
    plot_spectra_3d(all_WL, all_normalized, ['3D Normalized Waterfall Layout (', type, ')']);

    handle_csv_export(all_WL, all_normalized, all_time_scaled, filenames);
end

% --- Helper for Inline Conditions ---
function val = ifempty(cond, tV, fV), if cond, val = tV; else, val = fV; end; end

% --- Sub-function: CSV Handling ---
function handle_csv_export(wl, norm_data, time_data, filenames)
    answer = questdlg('Export processed data to two CSV files?', 'Export Data', 'Yes', 'No', 'No');
    if strcmp(answer, 'Yes')
        [save_name, save_path] = uiputfile('Spectra_MaxUnity.csv', 'Save Normalized Data As');
        if ~isequal(save_name, 0)
            export_csv_with_headers(fullfile(save_path, save_name), wl, norm_data, filenames);
            time_scaled_name = strrep(save_name, 'MaxUnity', 'TimeCorrected');
            if strcmp(time_scaled_name, save_name); time_scaled_name = ['TimeCorrected_', save_name]; end
            export_csv_with_headers(fullfile(save_path, time_scaled_name), wl, time_data, filenames);
        end
    end
end

function export_csv_with_headers(full_file_path, wl, data_matrix, headers)
    fid = fopen(full_file_path, 'w');
    fprintf(fid, 'Wavelength_nm');
    for i = 1:length(headers); fprintf(fid, ',%s', headers{i}); end
    fprintf(fid, '\n');
    fclose(fid);
    dlmwrite(full_file_path, [wl, data_matrix], '-append', 'delimiter', ',');
end

function plot_spectra(x, y, title_str, legs)
    figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
    plot(x, y, 'LineWidth', 2); 
    grid on;
    xlabel('Wavelength (nm)'); ylabel('Intensity');
    title(title_str);
    legend(legs, 'Interpreter', 'none', 'Location', 'eastoutside');
end

% --- New Sub-function: 3D Layout ---
function plot_spectra_3d(x, y, title_str)
    figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
    
    % Enumerate depth axis (file index tracking)
    z = 1:size(y, 2); 
    
    % waterfall needs matrices for coordinates or matching array layouts
    % Transposing y so rows map to depth (z) and columns map to wavelength (x)
    waterfall(x, z, y');
    
    % Enhancing aesthetics for line clarity in 3D
    view(45, 50); % Tilt perspective to easily view across depth layers
    grid on;
    xlabel('Wavelength (nm)', 'FontWeight', 'bold');
    ylabel('File Enumeration (Depth)', 'FontWeight', 'bold');
    zlabel('Intensity', 'FontWeight', 'bold');
    title(title_str, 'FontSize', 12);
    
    % Enhances visibility of individual line profiles over faces
    colormap([0 0.4470 0.7410]); % Clean default blue profile colors
end


% function batch_analyze_spectra_multi_type()
%     close all
%     % --- Setup Paths ---
%     default_dir = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_05_18__LabData\SHG Spectrum v Alignment (RunH SD66 id23.3r Post-Polish1)\TiSapphire\varyY';
%     default_dir = 'C:\Users\brent\Downloads\NichiaGC1_20C_450mA_4hr_hold05212026_VBG_spectra';
%     if ~exist(default_dir, 'dir'); default_dir = pwd; end
%     
%     % 1. Multi-select files
%     cd(default_dir);
%     [filenames, start_path] = uigetfile({'*.txt', 'Spectra Files (*.txt)'}, ...
%         'Select Spectra Files', default_dir, 'MultiSelect', 'on');
%     
%     if isequal(filenames, 0); return; end
%     if ischar(filenames); filenames = {filenames}; end 
%     
%     % --- CONFIGURATION ---
%     type = 'OceanOptics'; % Toggle between 'Avantes' or 'OceanOptics'
% %     type = 'Avantes'; % Toggle between 'Avantes' or 'OceanOptics'
%     
%     smooth_ans = questdlg('Apply low-pass smoothing (Savitzky-Golay)?', ...
%                          'Filter Settings', 'Yes', 'No', 'No');
%     do_smooth = strcmp(smooth_ans, 'Yes');
% 
%     num_files = length(filenames);
%     target_int_time = 100; % Target in ms
%     
%     % Initialize containers
%     all_WL = []; all_normalized = []; all_time_scaled = [];
%     labels = cell(1, num_files);
%     fit_center_wl = zeros(1, num_files);
%     centroid_wl = zeros(1, num_files); 
% 
%     fprintf('\n--- Processing %s Data ---\n', type);
%     
%     for k = 1:num_files
%         full_path = fullfile(start_path, filenames{k});
%         labels{k} = strrep(filenames{k}, '_', ' ');
% 
%         % --- Adaptive Integration Time Extraction ---
%         int_time = 1; 
%         fid = fopen(full_path, 'r');
%         search_limit = ifempty(strcmp(type, 'Avantes'), 9, 17);
%         for i = 1:search_limit
%             line_str = fgetl(fid);
%             if contains(lower(line_str), 'integration time')
%                 if strcmp(type, 'Avantes')
%                     nums = sscanf(line_str, '%*[^0-9.]%f');
%                     if ~isempty(nums); int_time = nums(1); end % Avantes usually in ms
%                 else % OceanOptics
%                     nums = sscanf(line_str, 'Integration Time (sec): %e');
%                     if ~isempty(nums); int_time = nums(1) * 1000; end % sec to ms
%                 end
%             end
%         end
%         fclose(fid);
% 
%         % --- Read Numeric Data (Your Updated Section) ---
%         switch type
%             case 'Avantes'
%                 data = dlmread(full_path, ';', 9, 0);
%             case 'OceanOptics'
%                 raw = importdata(full_path, '\t', 17);
%                 if isstruct(raw), data = raw.data; else, data = raw; end
%         end
%         
%         WL = data(:,1);
%         Intensity = data(:,2);
%         
%         if do_smooth
%             Intensity = sgolayfilt(Intensity, 3, 11);
%         end
% 
%         if k == 1; all_WL = WL; end 
% 
%         % --- Gaussian Fit (430-470nm) ---
%         idx = find(WL >= 430 & WL <= 470);
%         if ~isempty(idx) && length(idx) > 3
%             x = WL(idx); y = Intensity(idx);
%             y_fit = y; y_fit(y_fit <= 0) = min(y_fit(y_fit > 0)); 
%             logY = log(y_fit);
%             coeffs = polyfit(x, logY, 2);
%             fit_center_wl(k) = -coeffs(2) / (2 * coeffs(1));
%         else
%             fit_center_wl(k) = NaN;
%         end
% 
%         % --- Centroid Calculation ---
%         baseMask = (WL >= 200) & (WL <= 300);
%         avg_noise = ifempty(any(baseMask), mean(Intensity(baseMask)), 0);
%         auto_thresh = avg_noise + (10 * ifempty(any(baseMask), std(Intensity(baseMask)), 10));
% 
%         if ~isempty(idx)
%             wl_seg = WL(idx); int_seg = Intensity(idx);
%             [maxVal, ~] = max(int_seg);
%             if maxVal > auto_thresh
%                 cutOff = avg_noise + (maxVal - avg_noise) * 0.5;
%                 pMask = int_seg >= cutOff;
%                 centroid_wl(k) = sum(wl_seg(pMask) .* int_seg(pMask)) / sum(int_seg(pMask));
%             else
%                 centroid_wl(k) = NaN;
%             end
%         else
%             centroid_wl(k) = NaN;
%         end
%         
%         fprintf('%-40s | G: %-7.2f | C: %-7.2f\n', filenames{k}, fit_center_wl(k), centroid_wl(k));
% 
%         % Rescaling
%         max_val = max(Intensity); if max_val <= 0; max_val = 1; end 
%         all_normalized(:, k) = Intensity / max_val;
%         all_time_scaled(:, k) = Intensity * (target_int_time / int_time);
%     end
% 
%     % --- Visualization & Export ---
%     plot_spectra(all_WL, all_normalized, ['Normalized (', type, ')'], labels);
%     plot_spectra(all_WL, all_time_scaled, sprintf('Time-Corrected (%gms)', target_int_time), labels);
%     handle_csv_export(all_WL, all_normalized, all_time_scaled, filenames);
% end
% 
% % --- Helper for Inline Conditions ---
% function val = ifempty(cond, tV, fV), if cond, val = tV; else, val = fV; end; end
% 
% % (Sub-functions plot_spectra and handle_csv_export remain the same as previous)
% % --- Sub-function: CSV Handling ---
% function handle_csv_export(wl, norm_data, time_data, filenames)
%     answer = questdlg('Export processed data to two CSV files?', 'Export Data', 'Yes', 'No', 'No');
%     if strcmp(answer, 'Yes')
%         [save_name, save_path] = uiputfile('Spectra_MaxUnity.csv', 'Save Normalized Data As');
%         if ~isequal(save_name, 0)
%             export_csv_with_headers(fullfile(save_path, save_name), wl, norm_data, filenames);
%             time_scaled_name = strrep(save_name, 'MaxUnity', 'TimeCorrected');
%             if strcmp(time_scaled_name, save_name); time_scaled_name = ['TimeCorrected_', save_name]; end
%             export_csv_with_headers(fullfile(save_path, time_scaled_name), wl, time_data, filenames);
%         end
%     end
% end
% 
% function export_csv_with_headers(full_file_path, wl, data_matrix, headers)
%     fid = fopen(full_file_path, 'w');
%     fprintf(fid, 'Wavelength_nm');
%     for i = 1:length(headers); fprintf(fid, ',%s', headers{i}); end
%     fprintf(fid, '\n');
%     fclose(fid);
%     dlmwrite(full_file_path, [wl, data_matrix], '-append', 'delimiter', ',');
% end
% 
% function plot_spectra(x, y, title_str, legs)
%     figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
%     plot(x, y, 'LineWidth', 2); 
%     grid on; %xlim([445,458]);
%     xlabel('Wavelength (nm)'); ylabel('Intensity');
%     title(title_str);
%     legend(legs, 'Interpreter', 'none', 'Location', 'eastoutside');
% end

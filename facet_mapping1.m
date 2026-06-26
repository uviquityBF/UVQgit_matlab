%% Waveguide Chip Automated Alignment Tool (with Touch-Up Origin Calibration)
clc; clear; close all;

default_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_05_18__LabData\facet_mapping';
default_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_05_26_LabData\gl04_suruga\facet_mapping';
cd(default_path);

%% 1. File Selection via GUI Prompts
fprintf('Please select your files using the pop-up windows...\n');

% --- Select Layout CSV ---
[layout_name, layout_folder] = uigetfile('*.csv', 'STEP 1: Select your NOMINAL LAYOUT CSV file');
if isequal(layout_name, 0), error('File selection canceled. Script aborted.'); end
layout_file = fullfile(layout_folder, layout_name);
fprintf('-> Loaded Layout: %s\n', layout_name);

% --- Select Measured Alignment CSV ---
[meas_name, meas_folder] = uigetfile('*.csv', 'STEP 2: Select your TRANSPOSED MEASURED ALIGNMENT CSV file');
if isequal(meas_name, 0), error('File selection canceled. Script aborted.'); end
measured_file = fullfile(meas_folder, meas_name);
fprintf('-> Loaded Measurements: %s\n', meas_name);

% --- Optional Origin Calibration (Fibers Kissing Point) ---
has_origin_cal = questdlg('Do you have a fiber-tip contact origin calibration file?', ...
	'Origin Calibration', 'Yes', 'No', 'No');

origin_offsets = zeros(6,1); % Default to zero if not supplied
if strcmp(has_origin_cal, 'Yes')
    [cal_name, cal_folder] = uigetfile('*.csv', 'STEP 3: Select your TOUCH-UP ORIGIN CALIBRATION file (Same Format)');
    if ~isequal(cal_name, 0)
        fid_cal = fopen(fullfile(cal_folder, cal_name), 'r');
        cal_lines = textscan(fid_cal, '%s', 'Delimiter', '\n');
        fclose(fid_cal);
        cal_lines = cal_lines{1};
        
        if length(cal_lines) < 7
            error('The calibration file must contain exactly 7 horizontal rows.');
        end
        
        % Parse rows 2-7 to extract the numeric calibration coordinates
        for row = 2:7
            numeric_cells = strsplit(cal_lines{row}, ',');
            origin_offsets(row-1) = str2double(numeric_cells{1});
        end
        fprintf('-> Applied Fiber-Tip Origin Calibration.\n');
    end
end

output_file = fullfile(layout_folder, 'predicted_targets.csv');   

%% 2. Read Layout CSV (Skip 4 Headers)
try
    opts_layout = detectImportOptions(layout_file);
    if isprop(opts_layout, 'VariableNamingRule'), opts_layout.VariableNamingRule = 'preserve';
    elseif isprop(opts_layout, 'PreserveVariableNames'), opts_layout.PreserveVariableNames = true; end
    opts_layout.NumHeaderLines = 4; 
    layout_table = readtable(layout_file, opts_layout);
catch
    layout_table = readtable(layout_file, 'HeaderLines', 4, 'ReadVariableNames', false);
end

raw_identifiers = layout_table{:, 1};
if isnumeric(raw_identifiers)
    identifiers = strtrim(cellstr(num2str(raw_identifiers, '%.6g')));
else
    identifiers = strtrim(cellstr(raw_identifiers));
end
layout_input  = layout_table{:, 2:4};  
layout_output = layout_table{:, 5:7}; 

%% 3. Robust Parse of Transposed Measured Points CSV
fid = fopen(measured_file, 'r');
if fid == -1, error('Could not open the measurement file.'); end
raw_lines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
raw_lines = raw_lines{1};

meas_ids_raw = strsplit(raw_lines{1}, ',');
meas_ids_raw = strtrim(meas_ids_raw(~cellfun(@isempty, meas_ids_raw)));
num_measured = length(meas_ids_raw);

meas_data_matrix = zeros(6, num_measured);
for row = 2:7
    numeric_cells = strsplit(raw_lines{row}, ',');
    for col = 1:num_measured
        meas_data_matrix(row-1, col) = str2double(numeric_cells{col});
    end
end

% Extract raw measurements
meas_in_raw  = meas_data_matrix(1:3, :); 
meas_out_raw = meas_data_matrix(4:6, :);

% --- APPLY FIBER-TIP TOUCH ORIGIN CALIBRATION ---
meas_in_cal  = (meas_in_raw - origin_offsets(1:3))';
meas_out_cal = (meas_out_raw - origin_offsets(4:6))';

fprintf('Found %d alignment point(s) in the alignment file.\n', num_measured);

%% 4. Remap Layout Axes to Global Stage Frame
% layout x --> stage z | layout y --> stage x | layout z --> stage y
layout_in_stage_frame  = [layout_input(:,2),  layout_input(:,3),  layout_input(:,1)];
layout_out_stage_frame = [layout_output(:,2), layout_output(:,3), layout_output(:,1)];

nominal_length = mean(sqrt(sum((layout_out_stage_frame - layout_in_stage_frame).^2, 2)));

%% 5. Process & Mirror Measured Coordinates
corresponding_layout_in  = zeros(num_measured, 3);
corresponding_layout_out = zeros(num_measured, 3);
measured_in_global       = zeros(num_measured, 3);
measured_out_global      = zeros(num_measured, 3);

layout_lengths   = zeros(num_measured, 1);
physical_lengths = zeros(num_measured, 1);

for i = 1:num_measured
    id_str = meas_ids_raw{i};
    id_str = regexprep(id_str, '^id', ''); 
    
    idx = find(strcmp(identifiers, id_str));
    if isempty(idx)
        error('Measured ID "%s" not found in the Layout CSV file.', id_str);
    end
    
    corresponding_layout_in(i, :)  = layout_in_stage_frame(idx, :);
    corresponding_layout_out(i, :) = layout_out_stage_frame(idx, :);
    
    measured_in_global(i, :) = meas_in_cal(i, :);
    measured_out_global(i, :) = [-meas_out_cal(i,1), meas_out_cal(i,2), -meas_out_cal(i,3)];
    
    layout_lengths(i)   = norm(layout_out_stage_frame(idx, :) - layout_in_stage_frame(idx, :));
    physical_lengths(i) = norm(measured_out_global(i, :) - measured_in_global(i, :));
end

% %% 6. Correct Uniform Chip Length Error
% delta_length = mean(physical_lengths - layout_lengths);
% true_length = nominal_length + delta_length;
% 
% layout_out_stage_frame(:, 2) = layout_out_stage_frame(:, 2) + delta_length;
% corresponding_layout_out(:, 2) = corresponding_layout_out(:, 2) + delta_length;

%% 6. Correct Uniform Chip Length Error
% Calculate nominal length directly from the specific matched layout item 
% to guarantee no text columns or empty rows accidentally pollute the array math.
idx_first_match = find(strcmp(identifiers, regexprep(meas_ids_raw{1}, '^id', '')));
nominal_length_safe = norm(layout_out_stage_frame(idx_first_match, :) - layout_in_stage_frame(idx_first_match, :));

delta_length = mean(physical_lengths - layout_lengths);
true_length = nominal_length_safe + delta_length;

layout_out_stage_frame(:, 2) = layout_out_stage_frame(:, 2) + delta_length;
corresponding_layout_out(:, 2) = corresponding_layout_out(:, 2) + delta_length;


%% 7. Mathematical Matrix Alignment (Adaptive)
P_layout   = [corresponding_layout_in; corresponding_layout_out];
Q_measured = [measured_in_global; measured_out_global];

if num_measured == 1
    T = Q_measured(1, :) - P_layout(1, :);
    v_layout = (P_layout(2, :) - P_layout(1, :)) / norm(P_layout(2, :) - P_layout(1, :));
    v_meas   = (Q_measured(2, :) - Q_measured(1, :)) / norm(Q_measured(2, :) - Q_measured(1, :));
    
    rot_axis = cross(v_layout, v_meas);
    if norm(rot_axis) < 1e-6
        R = eye(3);
    else
        rot_axis = rot_axis / norm(rot_axis);
        rot_angle = acos(dot(v_layout, v_meas));
        K = [0, -rot_axis(3), rot_axis(2); rot_axis(3), 0, -rot_axis(1); -rot_axis(2), rot_axis(1), 0];
        R = eye(3) + sin(rot_angle)*K + (1-cos(rot_angle))*K^2;
    end
else
    centroid_P = mean(P_layout, 1);
    centroid_Q = mean(Q_measured, 1);
    
    P_centered = P_layout - centroid_P;
    Q_centered = Q_measured - centroid_Q;
    
    H = P_centered' * Q_centered;
    [U, ~, V] = svd(H);
    R = V * U';
    
    if det(R) < 0
        V(:, 3) = V(:, 3) * -1;
        R = V * U';
    end
    T = centroid_Q - (centroid_P * R');
end

%% 8. Map to Global Frame and Revert Output Stage Rotations & Origins
pred_global_in  = (layout_in_stage_frame * R') + T;
pred_global_out = (layout_out_stage_frame * R') + T;

pred_out_stage_cal = [-pred_global_out(:,1), pred_global_out(:,2), -pred_global_out(:,3)];

predicted_input_stage  = pred_global_in + origin_offsets(1:3)';
predicted_output_stage = pred_out_stage_cal + origin_offsets(4:6)';

%% 9. Format Data and Write to New CSV File
output_table = table(identifiers, ...
    predicted_input_stage(:,1),  predicted_input_stage(:,2),  predicted_input_stage(:,3), ...
    predicted_output_stage(:,1), predicted_output_stage(:,2), predicted_output_stage(:,3), ...
    'VariableNames', {'Identifier', 'Input_Stage_X', 'Input_Stage_Y', 'Input_Stage_Z', ...
                                    'Output_Stage_X', 'Output_Stage_Y', 'Output_Stage_Z'});
writetable(output_table, output_file);

%% 10. Extract Angles for Display
try
    euler_deg = rotm2eul(R) * (180/pi); 
catch
    pitch = -asin(R(3,1));
    yaw   = atan2(R(2,1)/cos(pitch), R(1,1)/cos(pitch));
    roll  = atan2(R(3,2)/cos(pitch), R(3,3)/cos(pitch));
    euler_deg = [yaw, pitch, roll] * (180/pi);
end

%% 11. 3D Scatter Plot Visualization
figure('Color', 'w', 'Position', [100, 100, 900, 650]);
hold on; grid on;

c_in_nom  = [0.2 0.5 0.8]; c_in_pred  = [0.1 0.3 0.6];
c_out_nom = [0.9 0.4 0.2]; c_out_pred = [0.7 0.2 0.1];

scatter3(layout_in_stage_frame(:,1), layout_in_stage_frame(:,2), layout_in_stage_frame(:,3), ...
    60, 'o', 'MarkeredgeColor', c_in_nom, 'LineWidth', 1.5, 'DisplayName', 'Nominal Inputs');
scatter3(pred_global_in(:,1), pred_global_in(:,2), pred_global_in(:,3), ...
    15, 'filled', 'MarkerFaceColor', c_in_pred, 'DisplayName', 'Predicted Inputs');
scatter3(layout_out_stage_frame(:,1), layout_out_stage_frame(:,2), layout_out_stage_frame(:,3), ...
    60, 'o', 'MarkeredgeColor', c_out_nom, 'LineWidth', 1.5, 'DisplayName', 'Nominal Outputs');
scatter3(pred_global_out(:,1), pred_global_out(:,2), pred_global_out(:,3), ...
    15, 'filled', 'MarkerFaceColor', c_out_pred, 'DisplayName', 'Predicted Outputs');

scatter3(measured_in_global(:,1), measured_in_global(:,2), measured_in_global(:,3), ...
    120, '^', 'filled', 'MarkerFaceColor', [0 0.6 0], 'MarkerEdgeColor', 'k', 'DisplayName', 'Measured Inputs');
scatter3(measured_out_global(:,1), measured_out_global(:,2), measured_out_global(:,3), ...
    120, 'v', 'filled', 'MarkerFaceColor', [0.6 0 0], 'MarkerEdgeColor', 'k', 'DisplayName', 'Measured Outputs');

axis equal; view(3); 
xlabel('Stage X (Calibrated Frame)'); ylabel('Stage Y (Calibrated Frame)'); zlabel('Stage Z (Calibrated Frame)');

title_str = {sprintf('True Chip Length: %.4f mm (Adjustment: %+.4f mm)', true_length, delta_length), ...
             sprintf('Calculated Tilts -> Yaw: %.3f° | Pitch: %.3f° | Roll: %.3f°', euler_deg(1), euler_deg(2), euler_deg(3))};
title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'bestoutside');
hold off;

%% Diagnostics Printout
fprintf('\n--- Processing Complete ---\n');
fprintf('Calculated Length Adjustment: %.4f mm\n', delta_length);
fprintf('Apparent True Length:         %.4f mm\n', true_length);
fprintf('Calculated System Tilts (Yaw, Pitch, Roll): [%.3f°, %.3f°, %.3f°]\n', euler_deg(1), euler_deg(2), euler_deg(3));
fprintf('Saved predicted coordinates to: %s\n\n', output_file);



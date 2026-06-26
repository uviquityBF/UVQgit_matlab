%% COMSOL Waveguide Mode Profile Visualizer
% Compatible with MATLAB R2018b and later

clear; clc; close all;
default_path = 'G:\Shared drives\uvq-Avo [external]\Gen 1 Product\NECSEL\mode_profiles_TM00_v_hornwidth';
cd(default_path);
%% COMSOL Waveguide Mode Profile Visualizer
% Robust Grid Interpolation Method (Compatible with all MATLAB versions)

clear; clc; close all;

%% 1. Select the COMSOL Text File
[filename, pathname] = uigetfile({'*.txt;*.dat', 'Text/Data Files (*.txt, *.dat)'}, ...
                                 'Select the COMSOL Export File');
if isequal(filename, 0) || isequal(pathname, 0)
    disp('User canceled the operation.');
    return;
end
fullpath = fullfile(pathname, filename);

%% 2. Read the Table Data, Skipping Comment Lines
opts = delimitedTextImportOptions('NumVariables', 4);
opts.Delimiter = ',';
opts.VariableTypes = {'double', 'double', 'double', 'double'};
opts.CommentStyle = '%';

dataTable = readtable(fullpath, opts);

rawX = dataTable.Var1;  % X coordinates (microns)
rawY = dataTable.Var2;  % Y coordinates (microns)
normE = dataTable.Var3; % Electric Field Norm (V/m)
powerZ = dataTable.Var4; % Power flow / Intensity (W/m^2)

% --- NEW: DATA CLEANING STEP ---
% Find rows where ANY coordinate or data value is NaN or Inf
badRows = isnan(rawX) | isinf(rawX) | ...
           isnan(rawY) | isinf(rawY) | ...
           isnan(normE) | isinf(normE) | ...
           isnan(powerZ) | isinf(powerZ);

% Keep only the clean data rows
rawX(badRows) = [];
rawY(badRows) = [];
normE(badRows) = [];
powerZ(badRows) = [];

% Inform the user if rows were removed
if any(badRows)
    fprintf('Data Cleaned: Removed %d rows containing NaN/Inf nodes from COMSOL export.\n', sum(badRows));
end
% -------------------------------

%% 3. Interpolate Scattered Data onto a Uniform 2D Grid
% Establish a clean target grid based on the min/max limits
gridPoints = 1024;
uniqueX = linspace(min(rawX), max(rawX), gridPoints);
uniqueY = linspace(min(rawY), max(rawY), gridPoints);

[X_grid, Y_grid] = meshgrid(uniqueX, uniqueY);

% Interpolate scattered data onto the clean mesh
% Now guaranteed to be free of Inf or NaN coordinates
E_matrix = griddata(rawX, rawY, normE, X_grid, Y_grid, 'linear');
P_matrix = griddata(rawX, rawY, powerZ, X_grid, Y_grid, 'linear');

% Re-fill any outside extrapolation edge values with zero
E_matrix(isnan(E_matrix)) = 0;
P_matrix(isnan(P_matrix)) = 0;
%% 4. Convert Spatial Units to mm (To match Zemax later)
x_axis_mm = uniqueX / 1000;
y_axis_mm = uniqueY / 1000;

%% 5. Plot the Waveguide Profiles
figure('Color', 'w', 'Position', [100, 100, 1200, 500]);

% Left Plot: Electric Field Norm
subplot(1, 2, 1);
imagesc(x_axis_mm, y_axis_mm, E_matrix);
set(gca, 'YDir', 'normal');
axis image; colormap(jet); colorbar;
grid on; set(gca, 'GridColor', 'w', 'GridAlpha', 0.3);
xlabel('X Position (mm)', 'FontSize', 11);
ylabel('Y Position (mm)', 'FontSize', 11);
title('Waveguide Mode: Electric Field Norm (V/m)', 'FontSize', 12);

% Right Plot: Power Flow / Intensity Profile
subplot(1, 2, 2);
imagesc(x_axis_mm, y_axis_mm, P_matrix);
set(gca, 'YDir', 'normal');
axis image; colormap(jet); colorbar;
grid on; set(gca, 'GridColor', 'w', 'GridAlpha', 0.3);
xlabel('X Position (mm)', 'FontSize', 11);
ylabel('Y Position (mm)', 'FontSize', 11);
title('Waveguide Mode: Power Flow Z-Component (W/m^2)', 'FontSize', 12);

% Main overall window title
annotation('textbox', [0, 0.92, 1, 0.08], ...
           'String', sprintf('COMSOL Export File: %s', filename), ...
           'EdgeColor', 'none', ...
           'HorizontalAlignment', 'center', ...
           'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
       

%% Master Waveguide Coupling & Mode Overlap Calculator
% Compatible with MATLAB R2018b and later

addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\waveguide');
clear; clc; close all;

% Setup default directory if it exists
default_path = 'G:\Shared drives\uvq-Avo [external]\Gen 1 Product\NECSEL';
if exist(default_path, 'dir'), cd(default_path); end

%% ==========================================
%% 1. LOAD AND PARSE ZEMAX POP DATA (LASER INPUT)
%% ==========================================
fprintf('--- Step 1: Loading Zemax POP Data ---\n');

% 1a. Load Zemax Irradiance File
[fileI, pathI] = uigetfile('*.txt', 'Select the Zemax POP IRRADIANCE File');
if isequal(fileI,0), disp('Operation canceled.'); return; end
POP_I = load_zemax_pop_txt(fullfile(pathI, fileI));

% 1b. Load Zemax Phase File
[fileP, pathP] = uigetfile('*.txt', 'Select the matching Zemax POP PHASE File');
if isequal(fileP,0), disp('Operation canceled.'); return; end
POP_P = load_zemax_pop_txt(fullfile(pathP, fileP));

% Grid Size & Spacing taken from the Irradiance file header (assumes Phase matches)
zGridX = POP_I.gridX; zGridY = POP_I.gridY;
I_zemax_raw = POP_I.data;
Phi_zemax_raw = POP_P.data;

% Zemax Axes (Centered at 0,0 in mm)
zX = POP_I.x_axis;
zY = POP_I.y_axis;
[Z_Xmesh, Z_Ymesh] = meshgrid(zX, zY);

%% ==========================================
%% 2. LOAD AND CLEAN COMSOL DATA (CHIP MODE)
%% ==========================================
fprintf('\n--- Step 2: Loading COMSOL Waveguide Data ---\n');

[fileC, pathC] = uigetfile({'*.txt;*.dat'}, 'Select the COMSOL Waveguide Export File');
if isequal(fileC,0), disp('Operation canceled.'); return; end

[cX_um, cY_um, cNormE, cPowerZ] = load_comsol_mode_profile(fullfile(pathC, fileC));
cX_mm = cX_um / 1000;   % Convert um to mm to align with Zemax
cY_mm = cY_um / 1000;   % Convert um to mm to align with Zemax

%% ==========================================
%% 3. DEFINE DYNAMIC DENSE INTEGRATION GRID (1% Threshold)
%% ==========================================
fprintf('\n--- Step 3: Determining 1%% Intensity Bounding Box ---\n');

% Find spatial limits where Zemax Intensity is > 1% of peak
maxI = max(I_zemax_raw(:));
maskI = I_zemax_raw >= (0.01 * maxI);
zX_active = Z_Xmesh(maskI);
zY_active = Z_Ymesh(maskI);

% Find spatial limits where COMSOL Power Intensity is > 1% of peak
maxP = max(cPowerZ);
maskC = cPowerZ >= (0.01 * maxP);
cX_active = cX_mm(maskC);
cY_active = cY_mm(maskC);

% Define the absolute outer envelope boundary encompassing both active domains
minX = min([min(zX_active), min(cX_active)]);
maxX = max([max(zX_active), max(cX_active)]);
minY = min([min(zY_active), min(cY_active)]);
maxY = max([max(zY_active), max(cY_active)]);

fprintf('Integration Domain Restricted To:\n');
fprintf('  X: [%.4f to %.4f] mm\n', minX, maxX);
fprintf('  Y: [%.4f to %.4f] mm\n', minY, maxY);

% Generate the pristine high-density uniform Integration Grid over this limited window
gridPoints = 512;
intX = linspace(minX, maxX, gridPoints);
intY = linspace(minY, maxY, gridPoints);
[Int_Xmesh, Int_Ymesh] = meshgrid(intX, intY);

%% ==========================================
%% 4. INTERPOLATION AND COMPLEX FIELD RECONSTRUCTION
%% ==========================================
fprintf('\n--- Step 4: Reconstructing Fields on Dense Integration Grid ---\n');

% Interpolate Zemax Raw Intensity and Phase onto the localized grid
I_int = griddata(zX, zY, I_zemax_raw, Int_Xmesh, Int_Ymesh, 'linear');
Phi_int = griddata(zX, zY, Phi_zemax_raw, Int_Xmesh, Int_Ymesh, 'linear');
I_int(isnan(I_int)) = 0; Phi_int(isnan(Phi_int)) = 0;

% Construct the Complex Electric Field Vector from the Laser Input
E_laser = sqrt(I_int) .* exp(1i * Phi_int);

% Interpolate COMSOL Scattered Mode Profile onto the exact same grid
E_chip = griddata(cX_mm, cY_mm, cNormE, Int_Xmesh, Int_Ymesh, 'linear');
E_chip(isnan(E_chip)) = 0;

%% ==========================================
%% 5. CALCULATE MODE COUPLING EFFICIENCY (OVERLAP INTEGRAL)
%% ==========================================
fprintf('\n--- Step 5: Calculating Overlap Integral ---\n');

% Numerator: |Double Integral( E_laser * conj(E_chip) dx dy) |^2
% Since grid spacing is uniform, dx*dy constants cancel out in the division.
numerator = abs(trapz(intY, trapz(intX, E_laser .* conj(E_chip), 2)))^2;

% Denominator: Integral(|E_laser|^2) * Integral(|E_chip|^2)
denom_laser = trapz(intY, trapz(intX, abs(E_laser).^2, 2));
denom_chip  = trapz(intY, trapz(intX, abs(E_chip).^2, 2));

coupling_efficiency = numerator / (denom_laser * denom_chip);

fprintf('\n======================================\n');
fprintf('  SUCCESS: OVERLAP COUPLING EFFICIENCY\n');
fprintf('  Calculated Coupling: %.2f%%\n', coupling_efficiency * 100);
fprintf('======================================\n');

%% ==========================================
%% 6. PLOT COMPARATIVE RESULTS
%% ==========================================
figure('Color', 'w', 'Position', [50, 50, 1400, 450]);

% Left Plot: Laser Field Magnitude on Integration Grid
subplot(1, 3, 1);
imagesc(intX*1000, intY*1000, abs(E_laser));
set(gca, 'YDir', 'normal'); axis image; colormap(jet); colorbar; grid on;
xlabel('X Position (�m)'); ylabel('Y Position (�m)');
title('Input Laser Field |E| (Interpolated)');

% Middle Plot: Waveguide Field Magnitude on Integration Grid
subplot(1, 3, 2);
imagesc(intX*1000, intY*1000, abs(E_chip));
set(gca, 'YDir', 'normal'); axis image; colormap(jet); colorbar; grid on;
xlabel('X Position (�m)'); ylabel('Y Position (�m)');
title('Waveguide Target Field |E| (Interpolated)');

% Right Plot: Laser Field Phase on Integration Grid
subplot(1, 3, 3);
imagesc(intX*1000, intY*1000, angle(E_laser));
set(gca, 'YDir', 'normal'); axis image; colormap(hsv); colorbar; grid on;
xlabel('X Position (�m)'); ylabel('Y Position (�m)');
title('Input Laser Phase Angle (rad)');

annotation('textbox', [0, 0.91, 1, 0.09], ...
           'String', sprintf('Master Overlap Efficiency Analysis\\nCalculated Coupling: %.2f%%', coupling_efficiency * 100), ...
           'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
           'FontSize', 14, 'FontWeight', 'bold');
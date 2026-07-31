%% modal_overlap.m
% Computes the power modal overlap integral between:
%   - A Zemax POP beam (irradiance + phase ASCII .txt files)
%   - A COMSOL waveguide/fiber mode profile (single-component CSV)
%
% Overlap definition (rigorous power overlap integral):
%
%         |  iint( E1 .* conj(E2) ) dA  |^2
%   eta = ----------------------------------
%         iint(|E1|^2 dA) * iint(|E2|^2 dA)
%
% where E1 is the complex Zemax field reconstructed from irradiance + phase,
% and E2 is the (real or complex) COMSOL mode field.
%
% Assumptions:
%   - Zemax irradiance file stores intensity I = |E|^2  =>  |E| = sqrt(I)
%   - Zemax phase file stores phase in radians
%   - COMSOL CSV has columns [x, y, field_value] with one header row
%   - Both fields are transverse; a single scalar component is used
%   - Both grids are interpolated onto a common Cartesian grid before integration
%
% Requires: MATLAB R2018b or later (no extra toolboxes needed)
%
% Author:  <your name>
% Date:    <date>

clear; clc; close all;

%% =========================================================
%  1.  USER SETTINGS  —  edit this block
%% =========================================================

% --- File paths ---
zemax_irradiance_file = 'irradiance.txt';   % Zemax POP irradiance ASCII file
zemax_phase_file      = 'phase.txt';        % Zemax POP phase ASCII file
comsol_csv_file       = 'mode_profile.csv'; % COMSOL field CSV  (x, y, F)

% --- Common interpolation grid ---
N_common   = 512;          % Number of points per side (NxN grid)
% Leave as [] to auto-set extents from the data overlap region,
% or specify manually as [x_min, x_max, y_min, y_max] in metres.
common_extent = [];        % e.g. [-50e-6, 50e-6, -50e-6, 50e-6]

% --- Zemax file format (number of header lines before data block) ---
% Standard Zemax POP ASCII export has a descriptive header followed by a
% row listing grid dimensions, then the data rows.  Adjust if your export
% differs.
zemax_header_lines = 16;   % lines to skip before numeric data

% --- Phase unwrapping ---
apply_phase_unwrap = true; % set false if phase is already continuous

% --- Plotting ---
show_plots = true;

%% =========================================================
%  2.  LOAD ZEMAX POP FILES
%% =========================================================
fprintf('Loading Zemax POP irradiance file: %s\n', zemax_irradiance_file);
[Iz, xz, yz] = load_zemax_ascii(zemax_irradiance_file, zemax_header_lines);

fprintf('Loading Zemax POP phase file:      %s\n', zemax_phase_file);
[Pz, ~,  ~]  = load_zemax_ascii(zemax_phase_file,      zemax_header_lines);

% Validate that irradiance and phase grids are consistent
assert(isequal(size(Iz), size(Pz)), ...
    'Irradiance and phase arrays must have the same dimensions.');

% Reconstruct complex field:  E = sqrt(I) * exp(i*phi)
if apply_phase_unwrap
    Pz = unwrap(unwrap(Pz, [], 1), [], 2);
end
E_zemax = sqrt(max(Iz, 0)) .* exp(1i .* Pz);

fprintf('  Zemax grid : %d x %d points,  x=[%.3g, %.3g] m,  y=[%.3g, %.3g] m\n', ...
    size(E_zemax,2), size(E_zemax,1), xz(1), xz(end), yz(1), yz(end));

%% =========================================================
%  3.  LOAD COMSOL CSV
%% =========================================================
fprintf('Loading COMSOL mode profile:       %s\n', comsol_csv_file);
[E_comsol_raw, xc, yc] = load_comsol_csv(comsol_csv_file);

fprintf('  COMSOL data: %d points,  x=[%.3g, %.3g] m,  y=[%.3g, %.3g] m\n', ...
    numel(xc), min(xc), max(xc), min(yc), max(yc));

%% =========================================================
%  4.  BUILD COMMON INTERPOLATION GRID
%% =========================================================
if isempty(common_extent)
    % Use the intersection (overlap) of both domains
    x_min = max(xz(1),   min(xc));
    x_max = min(xz(end), max(xc));
    y_min = max(yz(1),   min(yc));
    y_max = min(yz(end), max(yc));
    if x_min >= x_max || y_min >= y_max
        error('Zemax and COMSOL spatial extents do not overlap. Check units or common_extent setting.');
    end
else
    x_min = common_extent(1);  x_max = common_extent(2);
    y_min = common_extent(3);  y_max = common_extent(4);
end

x_common = linspace(x_min, x_max, N_common);
y_common = linspace(y_min, y_max, N_common);
[Xc, Yc] = meshgrid(x_common, y_common);

dx = x_common(2) - x_common(1);
dy = y_common(2) - y_common(1);

fprintf('\nCommon grid: %d x %d,  x=[%.3g, %.3g] m,  y=[%.3g, %.3g] m\n', ...
    N_common, N_common, x_min, x_max, y_min, y_max);

%% =========================================================
%  5.  INTERPOLATE BOTH FIELDS ONTO COMMON GRID
%% =========================================================

% --- Zemax (regular grid -> regular grid via interp2) ---
[Xz, Yz] = meshgrid(xz, yz);
E1_re = interp2(Xz, Yz, real(E_zemax), Xc, Yc, 'linear', 0);
E1_im = interp2(Xz, Yz, imag(E_zemax), Xc, Yc, 'linear', 0);
E1 = E1_re + 1i .* E1_im;

% --- COMSOL (scattered -> regular grid via griddata) ---
% COMSOL CSVs are often on a regular grid; griddata handles both cases.
E2 = griddata(xc, yc, E_comsol_raw, Xc, Yc, 'linear');
E2(isnan(E2)) = 0;    % outside convex hull -> zero (evanescent tails)

fprintf('Interpolation complete.\n');

%% =========================================================
%  6.  COMPUTE MODAL OVERLAP INTEGRAL
%% =========================================================
%
%  eta = | sum_ij( E1_ij * conj(E2_ij) ) * dx*dy |^2
%        -----------------------------------------------
%        sum_ij(|E1_ij|^2)*dx*dy  *  sum_ij(|E2_ij|^2)*dx*dy
%
%  The dx*dy area element appears in numerator squared and denominator
%  (each term once), so it cancels — but we keep it explicit for clarity
%  and dimensional correctness.

dA = dx * dy;

numerator   = abs( sum(sum( E1 .* conj(E2) )) * dA )^2;
norm_E1     = sum(sum( abs(E1).^2 )) * dA;
norm_E2     = sum(sum( abs(E2).^2 )) * dA;
denominator = norm_E1 * norm_E2;

if denominator == 0
    error('One or both fields have zero power on the common grid. Check units or spatial overlap.');
end

eta = numerator / denominator;

fprintf('\n========================================\n');
fprintf('  Modal Overlap  eta = %.6f  (%.4f %%)\n', eta, eta * 100);
fprintf('========================================\n');
fprintf('  Zemax field power  (on common grid): %.4g  [arb. units * m^2]\n', norm_E1);
fprintf('  COMSOL field power (on common grid): %.4g  [arb. units * m^2]\n', norm_E2);

%% =========================================================
%  7.  PLOTS
%% =========================================================
if show_plots
    x_um = x_common * 1e6;   % convert to microns for display
    y_um = y_common * 1e6;

    figure('Name', 'Modal Overlap — Field Comparison', 'NumberTitle', 'off', ...
           'Position', [100 100 1200 800]);

    subplot(2,3,1);
    imagesc(x_um, y_um, abs(E1).^2);
    axis image; colorbar; colormap(gca, 'hot');
    title('Zemax: Irradiance |E_1|^2');
    xlabel('x (\mum)'); ylabel('y (\mum)');

    subplot(2,3,2);
    imagesc(x_um, y_um, angle(E1));
    axis image; colorbar; colormap(gca, 'hsv');
    title('Zemax: Phase \angle E_1 (rad)');
    xlabel('x (\mum)'); ylabel('y (\mum)');

    subplot(2,3,3);
    imagesc(x_um, y_um, abs(E2).^2);
    axis image; colorbar; colormap(gca, 'hot');
    title('COMSOL: Mode Intensity |E_2|^2');
    xlabel('x (\mum)'); ylabel('y (\mum)');

    subplot(2,3,4);
    % Overlap integrand (amplitude)
    overlap_integrand = abs(E1 .* conj(E2));
    imagesc(x_um, y_um, overlap_integrand);
    axis image; colorbar; colormap(gca, 'parula');
    title('Overlap integrand |E_1 \cdot E_2^*|');
    xlabel('x (\mum)'); ylabel('y (\mum)');

    subplot(2,3,5);
    % 1-D cross-sections through centre
    mid = round(N_common/2);
    plot(x_um, abs(E1(mid,:)).^2 ./ max(abs(E1(mid,:)).^2 + eps), 'b-',  'LineWidth', 1.5); hold on;
    plot(x_um, abs(E2(mid,:)).^2 ./ max(abs(E2(mid,:)).^2 + eps), 'r--', 'LineWidth', 1.5);
    legend('Zemax E_1 (norm.)', 'COMSOL E_2 (norm.)', 'Location', 'best');
    xlabel('x (\mum)'); ylabel('Norm. Intensity');
    title('Horizontal cross-section (y = centre)');
    grid on;

    subplot(2,3,6);
    plot(y_um, abs(E1(:,mid)).^2 ./ max(abs(E1(:,mid)).^2 + eps), 'b-',  'LineWidth', 1.5); hold on;
    plot(y_um, abs(E2(:,mid)).^2 ./ max(abs(E2(:,mid)).^2 + eps), 'r--', 'LineWidth', 1.5);
    legend('Zemax E_1 (norm.)', 'COMSOL E_2 (norm.)', 'Location', 'best');
    xlabel('y (\mum)'); ylabel('Norm. Intensity');
    title('Vertical cross-section (x = centre)');
    grid on;

    sgtitle(sprintf('Modal Overlap  \\eta = %.4f  (%.2f %%)', eta, eta*100), ...
            'FontSize', 14, 'FontWeight', 'bold');
end

%% =========================================================
%  LOCAL FUNCTIONS
%% =========================================================

function [data, x_vec, y_vec] = load_zemax_ascii(filename, n_header)
% LOAD_ZEMAX_ASCII  Parse a Zemax POP Grid Sag / Irradiance / Phase ASCII file.
%
% Zemax POP ASCII layout (standard export):
%   Lines 1..n_header  : text header (title, wavelength, grid info, etc.)
%   Remaining lines    : Ny rows of Nx whitespace-delimited values
%
% The physical extent is read from the header line containing
% "X-grid spacing" and "Y-grid spacing" (metres).  If that pattern is not
% found the function falls back to unit pixel spacing with a warning.
%
% Returns:
%   data  [Ny x Nx] double array
%   x_vec [1 x Nx] coordinate vector (metres)
%   y_vec [Ny x 1] coordinate vector (metres)

    fid = fopen(filename, 'r');
    if fid < 0
        error('Cannot open file: %s', filename);
    end

    header_text = '';
    for k = 1:n_header
        line = fgetl(fid);
        if ischar(line)
            header_text = [header_text, line, newline]; %#ok<AGROW>
        end
    end

    % Read remaining numeric data
    raw = textscan(fid, '%f', 'CollectOutput', true);
    fclose(fid);
    raw = raw{1};

    % --- Parse grid size from header ---
    % Look for a line like:  "   64   64   ..." or "NX=  64, NY=  64"
    % Strategy: infer square root, or look for explicit Nx/Ny tokens.
    nx = []; ny = [];
    tokens = regexp(header_text, 'Grid Size\s*[:\=]?\s*(\d+)\s*[xX\,\s]+(\d+)', 'tokens', 'once');
    if ~isempty(tokens)
        nx = str2double(tokens{1});
        ny = str2double(tokens{2});
    else
        % Fallback: assume square grid
        n = round(sqrt(numel(raw)));
        if n^2 == numel(raw)
            nx = n; ny = n;
        else
            error(['Could not determine grid dimensions from header of %s.\n' ...
                   'Please adjust zemax_header_lines or extend the parser.'], filename);
        end
    end

    data = reshape(raw, nx, ny)';   % [ny x nx]

    % --- Parse grid spacing from header ---
    dx = []; dy = [];
    tok_dx = regexp(header_text, '[Xx][- _]?[Gg]rid\s*[Ss]pacing\s*[:\=]?\s*([\d\.eE\+\-]+)', 'tokens', 'once');
    tok_dy = regexp(header_text, '[Yy][- _]?[Gg]rid\s*[Ss]pacing\s*[:\=]?\s*([\d\.eE\+\-]+)', 'tokens', 'once');
    if ~isempty(tok_dx), dx = str2double(tok_dx{1}); end
    if ~isempty(tok_dy), dy = str2double(tok_dy{1}); end

    if isempty(dx) || isempty(dy)
        warning('Could not parse grid spacing from %s — using unit pixel spacing.\nCheck zemax_header_lines or extend the header parser.', filename);
        dx = 1; dy = 1;
    end

    x_vec = ((0:nx-1) - (nx-1)/2) * dx;   % centred coordinate
    y_vec = ((0:ny-1) - (ny-1)/2) * dy;
end


function [field, x, y] = load_comsol_csv(filename)
% LOAD_COMSOL_CSV  Load a COMSOL exported CSV with columns [x, y, field].
%
% Expected format:
%   Row 1  : header (column names — skipped)
%   Row 2+ : numeric data,  comma-separated
%
% Returns column vectors x, y and field values.
% The field may be real (magnitude) or complex; if COMSOL exports
% separate Re/Im columns, extend this function accordingly.

    fid = fopen(filename, 'r');
    if fid < 0
        error('Cannot open file: %s', filename);
    end
    fgetl(fid);   % discard header row
    raw = textscan(fid, '%f %f %f', 'Delimiter', ',', 'CollectOutput', true);
    fclose(fid);

    M = raw{1};
    if size(M, 2) < 3
        error('COMSOL CSV must have at least 3 columns: x, y, field_value.');
    end

    x     = M(:, 1);
    y     = M(:, 2);
    field = M(:, 3);   % extend here for complex: field = M(:,3) + 1i*M(:,4)
end
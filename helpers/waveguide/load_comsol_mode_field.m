function M = load_comsol_mode_field(filepath, unit_scale_to_um)
% Loads a COMSOL 6-column complex mode-field export:
%     x, y, Re(Ex), Im(Ex), Re(Ey), Im(Ey)
% comma delimited, comment/header lines starting with '%'. Rows containing
% NaN/Inf are dropped.
%
% USAGE:
%   M = load_comsol_mode_field(filepath);                % x,y already in um
%   M = load_comsol_mode_field(filepath, 1e-3);          % x,y exported in nm
%
% Returns a struct:
%   M.x, M.y   coordinates in um (column vectors)
%   M.Ex, M.Ey complex transverse field components (column vectors)
%
% Signed/complex fields are required (not NormE) so that sign changes in
% higher-order modes are preserved and orthogonality holds.

    if nargin < 2 || isempty(unit_scale_to_um), unit_scale_to_um = 1; end

    opts = delimitedTextImportOptions('NumVariables', 6);
    opts.Delimiter = ',';
    opts.VariableTypes = repmat({'double'}, 1, 6);
    opts.CommentStyle = '%';
    T = readtable(filepath, opts);

    D = T{:, 1:6};
    bad = any(isnan(D) | isinf(D), 2);
    D(bad, :) = [];
    if isempty(D)
        error(['No valid rows in %s. Expected 6 comma-separated columns: ' ...
               'x, y, Re(Ex), Im(Ex), Re(Ey), Im(Ey).'], filepath);
    end

    M.x  = D(:, 1) * unit_scale_to_um;
    M.y  = D(:, 2) * unit_scale_to_um;
    M.Ex = D(:, 3) + 1i * D(:, 4);
    M.Ey = D(:, 5) + 1i * D(:, 6);
end

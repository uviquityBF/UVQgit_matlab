function S = load_comsol_mode_set(filepath, unit_scale_to_um)
% Loads ONE COMSOL "Export > Data" spreadsheet (.csv) containing several mode
% solutions. Layout expected:
%     X, Y, then for each solution (mode), 4 columns:
%       real(<phys>.Ex), imag(<phys>.Ex), real(<phys>.Ey), imag(<phys>.Ey)
% Column headers must carry the solution tag after '@', e.g.
%       real(ewfd.Ex) (V/m) @ 2.1151 1     -> neff = 2.1151
% Header lines start with '%'. Rows containing NaN/Inf are dropped.
%
% USAGE:
%   S = load_comsol_mode_set(filepath);            % length unit read from header
%   S = load_comsol_mode_set(filepath, 1e-3);      % force units (nm -> um)
%
% Returns struct array S(1..nModes), in FILE order, with fields:
%   x, y        coordinates in um (column vectors, shared by all modes)
%   Ex, Ey      complex transverse fields (column vectors)
%   neff        real part of the effective index parsed from the header tag (NaN if not found)
%   neff_imag   imaginary part if the tag is complex (leaky mode), else 0
%   hasEy       false if the imag(Ey) column was missing/duplicated in the export
%               (then Ey is set to zero and must not be used)

    % ---- header ----
    fid = fopen(filepath, 'r');
    if fid == -1, error('Could not open %s', filepath); end
    nHeader = 0;  lastHdr = '';  unitLine = '';
    while true
        ln = fgetl(fid);
        if ~ischar(ln) || isempty(ln) || ln(1) ~= '%', break; end
        nHeader = nHeader + 1;  lastHdr = ln;
        if strncmpi(ln, '% Length unit', 13), unitLine = ln; end
    end
    fclose(fid);

    cols = strtrim(strsplit(lastHdr(2:end), ','));        % drop leading '%'
    nCols = numel(cols);
    if nCols < 6 || mod(nCols - 2, 4) ~= 0
        error(['Expected 2 + 4*nModes columns (X, Y, then Re/Im Ex, Re/Im Ey per mode); ' ...
               'header in %s has %d.'], filepath, nCols);
    end
    nM = (nCols - 2) / 4;

    if nargin < 2 || isempty(unit_scale_to_um)
        u = '';  tk = regexp(unitLine, 'Length unit,\s*(.*)$', 'tokens', 'once');
        if ~isempty(tk), u = strtrim(tk{1}); end
        if strcmp(u, 'nm'), unit_scale_to_um = 1e-3;
        elseif strcmp(u, 'mm'), unit_scale_to_um = 1e3;
        elseif strcmp(u, 'm'), unit_scale_to_um = 1e6;
        else, unit_scale_to_um = 1; end                   % um (or unreadable micro sign)
    end

    % ---- data ----
    D = dlmread(filepath, ',', nHeader, 0);
    D(any(isnan(D) | isinf(D), 2), :) = [];
    if isempty(D) || size(D, 2) ~= nCols
        error('Data block in %s has %d columns, header lists %d.', filepath, size(D, 2), nCols);
    end

    S = struct('x', cell(nM, 1), 'y', [], 'Ex', [], 'Ey', [], 'neff', [], 'neff_imag', [], 'hasEy', []);
    x = D(:, 1) * unit_scale_to_um;
    y = D(:, 2) * unit_scale_to_um;
    for k = 1:nM
        c = 2 + 4 * (k - 1);                              % column offset for this mode
        h = cols(c + 1:c + 4);
        if isempty(regexp(h{1}, '^real\(.*\.Ex\)', 'once')) || ...
           isempty(regexp(h{2}, '^imag\(.*\.Ex\)', 'once'))
            error('Mode %d: columns 1-2 should be real/imag(...Ex), found "%s", "%s".', k, h{1}, h{2});
        end
        hasEy = ~isempty(regexp(h{3}, '^real\(.*\.Ey\)', 'once')) && ...
                ~isempty(regexp(h{4}, '^imag\(.*\.Ey\)', 'once'));
        % solution tag: "2.1151" or complex "1.7773-3.2622E-4i" (leaky mode)
        tk = regexp(h{1}, '@\s*([0-9.]+)(?:([+-][0-9.]+(?:[Ee][+-]?\d+)?)i)?', 'tokens', 'once');

        S(k).x = x;  S(k).y = y;
        S(k).Ex = D(:, c + 1) + 1i * D(:, c + 2);
        if hasEy
            S(k).Ey = D(:, c + 3) + 1i * D(:, c + 4);
        else
            S(k).Ey = zeros(size(x));
        end
        S(k).hasEy = hasEy;
        if isempty(tk), S(k).neff = NaN; else, S(k).neff = str2double(tk{1}); end
        S(k).neff_imag = 0;
        if ~isempty(tk) && numel(tk) > 1 && ~isempty(tk{2}), S(k).neff_imag = str2double(tk{2}); end
    end
    if any(~[S.hasEy])
        warning('load_comsol_mode_set:noEy', ...
            ['%d of %d modes lack imag(Ey) in the export (4th expression per mode is not ' ...
             'imag(...Ey)). Ey is unusable; re-export with imag(ewfd.Ey).'], sum(~[S.hasEy]), nM);
    end
end

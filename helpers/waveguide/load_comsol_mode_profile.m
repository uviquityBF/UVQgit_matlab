function [x, y, normE, powerZ] = load_comsol_mode_profile(filepath)
% Loads a COMSOL 4-column mode-profile export (x, y, NormE, PowerZ), comma
% delimited with comment-style header lines starting with '%'.
% Rows containing NaN/Inf in any column are dropped.
%
% USAGE:
%   [x, y, normE, powerZ] = load_comsol_mode_profile(filepath);
%
% x, y are in whatever units the COMSOL export used (typically microns) --
% unit conversion is left to the caller.

    opts = delimitedTextImportOptions('NumVariables', 4);
    opts.Delimiter = ',';
    opts.VariableTypes = {'double', 'double', 'double', 'double'};
    opts.CommentStyle = '%';

    dataTable = readtable(filepath, opts);

    x = dataTable.Var1;
    y = dataTable.Var2;
    normE = dataTable.Var3;
    powerZ = dataTable.Var4;

    badRows = isnan(x) | isinf(x) | isnan(y) | isinf(y) | ...
              isnan(normE) | isinf(normE) | isnan(powerZ) | isinf(powerZ);
    x(badRows) = [];
    y(badRows) = [];
    normE(badRows) = [];
    powerZ(badRows) = [];

    if any(badRows)
        fprintf('Data Cleaned: Removed %d rows containing NaN/Inf nodes from COMSOL export (%s).\n', ...
            sum(badRows), filepath);
    end
end

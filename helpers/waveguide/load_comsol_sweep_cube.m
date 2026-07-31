function S = load_comsol_sweep_cube(filepath, sheet, range, varargin)
% Loads a COMSOL mode-index parameter-sweep table (xlsread) and reshapes it
% into a 3D [mode x P1 x P2] cube. Two input table layouts are supported
% via the 'Format' name-value option:
%
%   'header_rows' (default) -- row 1 of the range holds the Param1 value
%     (e.g. wavelength) for every column, row 2 holds the Param2 value
%     (e.g. width) for every column, rows 3:end are modes. Columns are
%     grouped by matching Param2 value.
%       S = load_comsol_sweep_cube(filepath, sheet, range)
%       S.cube          [Nmodes x Nvals1 x Nvals2]
%       S.P1s, S.P2s    raw header rows
%       S.P1_u, S.P2_u  unique Param1 / Param2 values
%       S.Nvals1, S.Nvals2, S.Nm
%
%   'column_groups' -- columns come pre-grouped into fixed-size blocks of
%     GroupSize columns, one block per Param2 value; row 1 within a block
%     holds the (scalar) Param2 value for the whole block, row 2 holds the
%     Param1 value for each column in the block, rows 3:end are modes.
%       S = load_comsol_sweep_cube(filepath, sheet, range, 'Format', 'column_groups', 'GroupSize', grp_n_col)
%       S.A             the raw xlsread matrix
%       S.Ngrps         number of groups found
%       S.M             1xNgrps cell array of raw group blocks
%       S.P2 (w)        Ngrps x 1 vector, one Param2 value per group
%       S.P1 (X)        1xNgrps cell array, Param1 vector per group (column vector)
%       S.Ys            1xNgrps cell array, [Param1 x mode] matrix per group

    p = inputParser;
    addParameter(p, 'Format', 'header_rows');
    addParameter(p, 'GroupSize', []);
    parse(p, varargin{:});
    fmt = p.Results.Format;
    grp_n_col = p.Results.GroupSize;

    A = xlsread(filepath, sheet, range);

    switch fmt
        case 'header_rows'
            P1s = A(1,:);
            P2s = A(2,:);
            P1s_u = unique(P1s(~isnan(P1s)));
            P2s_u = unique(P2s(~isnan(P2s)));
            Nvals1 = length(P1s_u);
            Nvals2 = length(P2s_u);
            Nm = size(A,1) - 2;

            cube = zeros(Nm, Nvals1, Nvals2);
            for k = 1:Nvals2
                p2 = P2s_u(k);
                icols_sel = find(P2s == p2);
                cube(:,:,k) = A(3:end, icols_sel);
            end

            S.cube = cube;
            S.P1s = P1s;
            S.P2s = P2s;
            S.P1_u = P1s_u;
            S.P2_u = P2s_u;
            S.Nvals1 = Nvals1;
            S.Nvals2 = Nvals2;
            S.Nm = Nm;

        case 'column_groups'
            if isempty(grp_n_col)
                error('load_comsol_sweep_cube: ''column_groups'' format requires the ''GroupSize'' option.');
            end
            Ngrps = floor(size(A,2) / grp_n_col);
            cstart = mod(size(A,2), grp_n_col);

            M = cell(1, Ngrps);
            w = zeros(Ngrps, 1);
            X = cell(1, Ngrps);
            Ys = cell(1, Ngrps);
            for kg = 1:Ngrps
                M{kg} = A(:, cstart + (grp_n_col*(kg-1)+1 : grp_n_col*kg));
                w(kg,1) = M{kg}(1,1);
                X{kg} = M{kg}(2,:)';
                Ys{kg} = M{kg}(3:end,:)';
            end

            S.A = A;
            S.M = M;
            S.P2 = w;
            S.P1 = X;
            S.Ys = Ys;
            S.Ngrps = Ngrps;

        otherwise
            error('load_comsol_sweep_cube: unknown Format ''%s''.', fmt);
    end
end

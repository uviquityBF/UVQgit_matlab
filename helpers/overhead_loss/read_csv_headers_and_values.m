function C = read_csv_headers_and_values(filepath)
% Reads any CSV file into a full string matrix via readtable/detectImportOptions.
%
%   C = read_csv_headers_and_values(filepath)
%
% C - string(table2cell(readtable(filepath))), i.e. every cell as a string,
%     in the same row/column layout as the source file (variable names
%     consumed by readtable are NOT included as a row of C).
%
% This is the boilerplate previously duplicated (identically) inside
% OverheadLoss_figure_Maker.m and overheadLoss_open_CSVs_and_compile3.m.
% The two callers use different CSV layouts downstream -- figure_Maker.m
% treats C as a wide table (one column per field, e.g. C(:,1)=Run,
% C(:,2)=wafer, ...), while overheadLoss_open_CSVs_and_compile3.m reads
% "*Params*.csv" files laid out as [~, headerName, value] triples and pulls
% headers = C(:,2); values = C(:,3) at the call site. Only this common
% read-and-stringify step is shared; column semantics stay with each caller.

    opts = detectImportOptions(filepath);
    T = readtable(filepath, opts);
    C = string(table2cell(T));
end

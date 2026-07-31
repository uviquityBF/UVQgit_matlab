function [fileStruct, filePaths, matData] = load_intermediate_mats(rootFolder)
% Recursively finds and loads all *intermediateData*.mat files under rootFolder.
%
%   [fileStruct, filePaths, matData] = load_intermediate_mats(rootFolder)
%
% fileStruct - dir() listing for every matched .mat file (folder, name, ...)
% filePaths  - cell array of full paths, fullfile(fileStruct(k).folder, fileStruct(k).name)
% matData    - cell array, one struct per file, holding every variable stored
%              in that .mat file (x_index_array, y_data_array, hdrNormalized, ...)
%              Kept as a cell array (not a struct array) since not every
%              file is guaranteed to contain identical variable sets.
%
% Shared batch-loading loop, factored out of OverheadLoss_IntermediateData.m
% and OverheadLoss_IntermediateData_DefectDetection.m.

    fileStruct = dir(fullfile(rootFolder, '**', '*.mat'));
    Nf = numel(fileStruct);
    fprintf('Found %d .mat files under %s\n', Nf, rootFolder);

    filePaths = cell(Nf,1);
    matData   = cell(Nf,1);
    for k = 1:Nf
        filePaths{k} = fullfile(fileStruct(k).folder, fileStruct(k).name);
        matData{k}   = load(filePaths{k});
        fprintf('  Loaded %d/%d: %s\n', k, Nf, fileStruct(k).name);
    end
end

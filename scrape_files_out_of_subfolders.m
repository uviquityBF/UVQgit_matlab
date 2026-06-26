%% File Extractor Script
% This script searches recursively for files matching a regex and 
% copies them to a 'copies' directory at the top level.


function scrape_files_out_of_subfolders()
    % --- Configuration ---
%     rootPath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Overhead Loss Historical Data Survey\Run C';           % The starting directory
%    rootPath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Overhead Loss Historical Data Survey\GaugeLots';
%     rootPath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run I\M3266';
    rootPath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\GaugeLots\GLI\M3266\SD01\Standardized Testing\Overhead Loss';
    pattern  = '_Overhead_Loss_Results - All';       % The regexp (e.g., matches all .txt files)
    destDir  = fullfile(rootPath, 'copies (all)');

    % --- Step 1: Create the destination folder ---
    if ~exist(destDir, 'dir')
        mkdir(destDir);
        fprintf('Created directory: %s\n', destDir);
    end

    % --- Step 2: Search for all files recursively ---
    % The '**' wildcard tells dir to look through all subfolders
    allFiles = dir(fullfile(rootPath, '**', '*.*'));

    % Remove directories from the list, keep only files
    allFiles = allFiles(~[allFiles.isdir]);

    % --- Step 3: Filter and Copy ---
    count = 0;
    for i = 1:length(allFiles)
        fileName = allFiles(i).name;
        sourcePath = fullfile(allFiles(i).folder, fileName);

        % Skip files already inside the 'copies' folder to avoid infinite loops
        if contains(allFiles(i).folder, destDir)
            continue;
        end

        % Check if filename matches the regular expression
        if ~isempty(regexp(fileName, pattern, 'once'))
            targetPath = fullfile(destDir, fileName);

            % Handle potential filename collisions (same name in different folders)
            if exist(targetPath, 'file')
                [~, name, ext] = fileparts(fileName);
                targetPath = fullfile(destDir, sprintf('%s_%d%s', name, count, ext));
            end

            [status, msg] = copyfile(sourcePath, targetPath);

            if status
                fprintf('Copied: %s\n', fileName);
                count = count + 1;
            else
                warning('Failed to copy %s: %s', fileName, msg);
            end
        end
    end

    fprintf('\nDone! Successfully copied %d files to %s\n', count, destDir);
    
end
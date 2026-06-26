%% scour all folder paths and pull images

function  compile_overhead_images()
    basePath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run C\AM123';
%     basePath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run F\AM136';
    basePath = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run F\AM136\SD43\R&R Study';
    pattern = '*_Overhead_Loss_Results - All.png';  % You can use wildcards:   '*report*.pdf'; 
    
    matching_filepaths = findMatchingFiles(basePath, pattern);

    disp('Matching files:');
    disp(matching_filepaths);
    
    %% open "ALL Results" and extract images
    if ~isdir([basePath,'\image_compilation\'])
        mkdir([basePath,'\image_compilation\']);
    end
    
    Nfiles = length(matching_filepaths);
    for kf=1:Nfiles
        disp(kf); disp(matching_filepaths{kf});
        I0 = imread(matching_filepaths{kf});
        
        Isub1 = I0(280:530,360:660);
        Isub2 = I0(6:30,850:1150);
        Inew = [Isub2;Isub1];        
        
        strpieces = split(matching_filepaths{kf},{'\','_Overhead_Loss_Results - All.png'});
        new_filename = [strpieces{8},'_',strpieces{9},'_','_',strpieces{12},'_',strpieces{13},'.png'];
        
        imwrite(Inew,[basePath,'\image_compilation\',new_filename]);
        
    end
   
    %% open "FIT Results" and extract Fit CUrves
    pattern2 = '*_Overhead_Loss_Results.png';  % You can use wildcards:   '*report*.pdf'; 
    
    matching_filepaths2 = findMatchingFiles(basePath, pattern2);
    
    if ~isdir([basePath,'\results_compilation\'])
        mkdir([basePath,'\results_compilation\']);
    end
    
    Nfiles = length(matching_filepaths2);
    for kf=1:Nfiles
        disp(kf); disp(matching_filepaths2{kf});
        I0 = imread(matching_filepaths2{kf});
        
        strpieces = split(matching_filepaths2{kf},{'\','_Overhead_Loss_Results.png'});
        new_filename = [strpieces{8},'_',strpieces{9},'_','_',strpieces{12},'_',strpieces{13},'.png'];
        
        imwrite(I0,[basePath,'\results_compilation\',new_filename]);
        
    end
    
end




function matchingFiles = findMatchingFiles(baseDir, matchPattern)
% FINDMATCHINGFILES Recursively search for files matching a string/pattern.
%
% matchingFiles = findMatchingFiles(baseDir, matchPattern)
%
% Inputs:
%   baseDir      - Root directory to search from
%   matchPattern - Pattern to match filenames (e.g., '*.txt', '*report*')
%
% Output:
%   matchingFiles - Cell array of full file paths that match the pattern

    % Initialize output
    matchingFiles = {};

    % Get list of all items in the current directory
    dirData = dir(baseDir);

    % Remove '.' and '..'
    dirData = dirData(~ismember({dirData.name}, {'.', '..'}));

    for i = 1:length(dirData)
        fullPath = fullfile(baseDir, dirData(i).name);

        if dirData(i).isdir
            % Recurse into subdirectories
            matchingFiles = [matchingFiles; findMatchingFiles(fullPath, matchPattern)];
        else
            % Match current file name to the pattern
            matchList = dir(fullfile(baseDir, matchPattern));
            matchNames = {matchList.name};

            if any(strcmp(dirData(i).name, matchNames))
                matchingFiles{end+1,1} = fullPath;
            end
        end
    end
end

% function matchingFiles = findMatchingFiles(baseDir, matchPattern)
% % FINDMATCHINGFILES Recursively search for files matching a string/pattern.
% %
% % matchingFiles = findMatchingFiles(baseDir, matchPattern)
% %
% % Inputs:
% %   baseDir      - Root directory to search from
% %   matchPattern - Pattern to match filenames (e.g., '*.txt', '*report*')
% %
% % Output:
% %   matchingFiles - Cell array of full file paths that match the pattern
% 
%     % Initialize output
%     matchingFiles = {};
% 
%     % Get list of all files and folders in the current directory
%     dirData = dir(baseDir);
% 
%     % Filter out . and .. directories
%     dirData = dirData(~ismember({dirData.name}, {'.', '..'}));
% 
%     for i = 1:length(dirData)
%         fullPath = fullfile(baseDir, dirData(i).name);
%         if dirData(i).isdir
%             % Recurse into subdirectory
%             matchingFiles = [matchingFiles; findMatchingFiles(fullPath, matchPattern)];
%         else
%             % Check if file name matches the pattern
%             
%             if ~isempty(dir(fullfile(baseDir, matchPattern))) && ...
%                     any( strcmp(dirData(i).name, {dir(fullfile(baseDir, matchPattern)).name}))
%                 matchingFiles{end+1,1} = fullPath;
%             end
%         end
%     end
% end
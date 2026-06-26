function rotate_images_batch
    dirs.root = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2026_02_09_LabData\GL01_id2.7 brightfield inspection';
    rootFolder = uigetdir(dirs.root, 'Select folder containing MAT files');
    fileStruct = dir(fullfile(rootFolder, '**', '*.png'));
    Nf = length(fileStruct)
    
    for kf = 1:Nf

        path_and_file = fullfile(fileStruct(kf).folder, fileStruct(kf).name);
        %split path and file
        s = split(path_and_file,'\');
        filename = s{end};
        tmp = split(path_and_file,filename);
        path = tmp{1};

        I0 = imread(path_and_file);
        
        I1 = imrotate(I0,180);
        
        imwrite(I1,[path,'/rotated/',filename]);
        
    end
   
end
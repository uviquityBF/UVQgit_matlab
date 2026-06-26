
% 1. Open a bunch of TIFF files
% 2. Carry out analysis by hand (draw from each section)

function TIFF_Analysis_SHG_from_WG()
%% get data
    close all;
    clear;
    dirs.f  = 'C:\Users\brent\MATLAB\UVQ';
    dirs.f = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    start_path = 'C:\Users\brent\MATLAB\matDat\';
    filepath = 'C:\Users\brent\MATLAB\matDat\';
    start_path  = 'C:\Users\brent\Pictures\';
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_11_03_LabData\UV_Camera -- SHG from WG';
    %Select all files within SubFolders
    %     toppath = uigetdir(start_path);
    toppath = 0;
    %Select by Hand
    if toppath==0
        cd(start_path);
        [filenames,start_path] = uigetfile({'*.tiff';'*.png'},'Multiselect','on');
        
         Nf = length(filenames);
         if isstr(filenames)
             filename = filenames; clear filenames;
             filenames{1} = filename;
             Nf=1;
         end
        for kf=1:Nf
            paths{kf} = start_path;
        end
    else
        cd(toppath);
%         A = dir('**/*.tiff');  %Search within sub folders
        A = dir('**/*.PNG');  %Search within sub folders
        Nf = length(A)
        for kf=1:Nf
            filenames{kf} = A(kf).name;
            paths{kf} = [A(kf).folder,'/'];
        end
    end
    
    
%% Do analyses in sections below ..

 %% [1] Ratio of Two Images
    for kf=1:Nf
        I0{kf} = imread([paths{kf},filenames{kf}]);
        if size(I0{kf},3)>1
            I0{kf} = rgb2gray(I0{kf}(:,:,1:3));
        end
%        figure; imshow(I0{kf});     
    end
    figure; imagesc(I0{kf});
    thresh=110;
    isel = find((I0{1}(:) > thresh)& (I0{2}(:) > thresh) & (I0{1}(:)<255));
    isel_not = setdiff( [1:size(I0{1},1)*size(I0{1},2)], isel);
    Iratio = double(I0{2})./double(I0{1});
    Iratio(isel_not) = 0;
    figure; subplot(1,3,1);  imshow(I0{1}); xlabel(replace(filenames{1},'_',' '));
    subplot(1,3,2);  imshow(I0{2});  xlabel(replace(filenames{2},'_',' '));
    subplot(1,3,3);  imagesc(Iratio); caxis([0.8,1.2]); colorbar;
    [yh,xh] = hist(Iratio(isel),[0.0:0.02:2]);
     figure; bar(xh,yh); xlabel('Ratio of Pixel Value (image 2 /  image 1) = 2xCT450rb / 1xCT450rb'); ylabel('number of pixels')
    title(['median = ',num2str(median(Iratio(isel)))]);
    
 %% Observe Spacing of other points
 kf = 1;
 figure;
 for kf=1:min(Nf,6)
    subplot(2,3,kf); imagesc(log10(double(I0{kf}-median(I0{kf}))));colorbar
 end
 
 
 %% Quantify Bright spot on image   
 % Loop Over Files
    for kf=1:Nf
        I0 = imread([paths{kf},filenames{kf}]);
        if size(I0,3)>1
            I0 = rgb2gray(I0(:,:,1:3));
        end
        % imshow(I0);
        

        %-------- TOTAL COUNTS ------------
        TotalImageCounts(kf,1) = sum(I0(:));
        BG_sample_region = I0(1:min([20,size(I0,1)]),1:min([20,size(I0,2)]));
        BG_imagecounts_estimate(kf,1)   = size(I0,1)*size(I0,2)* mean(mean(BG_sample_region));
        BG_stdev_estimate(kf,1)   =  std(double(BG_sample_region(:)));
        TOT_stdev_estimate(kf,1)   =  std(double(I0(:)));

        Nrows =  size(I0,1);  Ncols = size(I0,2);
        %-------- THRESHOLD IMAGE & CALCULATE COUNTS above BG------------


        %find Peak (max value)
        [maxval(kf),imax] = max(I0(:));
        [rmax,cmax] = ind2sub(size(I0),imax);

        %ROI IMAGE -- around the peak
%          I0roi = I0(rmax-200:rmax+200,cmax-200:cmax+200);
         I0roi = I0;


        scale=2.0;  %if(kf>2) scale = 1.5;   else scale=0.9; end;
        [level,EM] = graythresh(I0roi);
%        I1 = imbinarize(I0roi,level*scale);
%        thresh = graythresh(I0roi);
        I1 = im2bw(I0roi,level*scale);

        %         threshlevels = multithresh(I0roi,4);
        %         I1 = imquantize(I0,threshlevels);
        %         figure; imshow(I1);

        isel = find(I1(:)==1);
        SumCounts(kf,1) = sum(I0roi(isel));
        NumPix(kf,1) = length(isel);


        %SHOW FIGURE
        [rsel,csel] = ind2sub(size(I0roi),isel);

%         if kf==1 figure; subplot(1,Nf,1)
%         else subplot(1,Nf,kf)
%         end
        warning off

        hf1 = figure;

       	imshow(I0roi); hold all; plot(csel,rsel,'.','Color','r');
        title(replace(filenames{kf},'_',' '));
%         saveas(hf1,[toppath,'/',replace(filenames{kf},'.tiff','selectedPoints.jpg')], 'jpg');

        %PRINT RESULT
        X = sprintf('%s , %d , %d, %d, %d, %d, %d, %d, %d, %d',...
            filenames{kf},SumCounts(kf,1),NumPix(kf,1),maxval(kf),...
            TotalImageCounts(kf,1),BG_imagecounts_estimate(kf,1), BG_stdev_estimate(kf,1),TOT_stdev_estimate(kf,1), Nrows, Ncols);
        disp(X);

        warning onTot


    end

     disp(maxval)




end


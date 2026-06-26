
% 1. Open TIFF containing a single peak
% 2. Threshold the image
% 3. Computue total counts in the spot around a single peak
% 4. Report Peak value (need to flag Saturated Pixels = 255_

function TIFF_Analysis_Counts_from_DOUBLE_Peak()
%% get data
    close all;
    clear;
    dirs.f  = 'C:\Users\brent\MATLAB\UVQ';
    dirs.f = 'G:\My Drive\Analyses(BF)\Matlab\UVQ';
    start_path = 'C:\Users\brent\MATLAB\matDat\';
    filepath = 'C:\Users\brent\MATLAB\matDat\';
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run C\AM123\SD30\directional couplers\rear';
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

 %define ROIs to divide Image
    ROIs{1}   = [400,50;700,500];    %LEFT,  [row,col; row,col]
    ROIs{2}   = [400,500;700,950];   %RIGHT, [row,col; row,col]
    
 % Loop Over Files
    for kf=1:Nf
        
        %load image
        I0 = imread([paths{kf},filenames{kf}]);
        if size(I0,3)>1
            I0 = rgb2gray(I0(:,:,1:3));
        end
        
        %show ROIs
        hf0 = figure; imshow(I0);
        x1=ROIs{1}(1,2); x2=ROIs{1}(2,2); y1=ROIs{1}(1,1); y2=ROIs{1}(2,1);
        hold all;  plot([x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1],'g');
        x1=ROIs{2}(1,2); x2=ROIs{2}(2,2); y1=ROIs{2}(1,1); y2=ROIs{2}(2,1);
        hold all;  plot([x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1],'g');
        pause;
        close(hf0);

        %-------- TOTAL COUNTS ------------
        BG_sample_region                = I0(1:min([20,size(I0,1)]),1:min([20,size(I0,2)]));
        BG_avg                          = mean(mean(BG_sample_region));
        BG_imagecounts_estimate(kf,1)   = size(I0,1)*size(I0,2)*BG_avg ;
        BG_stdev_estimate(kf,1)         =  std(double(BG_sample_region(:)));
        TotalImageCounts(kf,1)          = sum(I0(:));
        TOT_stdev_estimate(kf,1)        =  std(double(I0(:)));
        maxval_fullimg                  = max(I0(:));

        Nrows =  size(I0,1);  Ncols = size(I0,2);
        %-------- THRESHOLD IMAGE & CALCULATE COUNTS above BG------------

        Nroi = length(ROIs);
        for kroi = 1:Nroi
            
            %ROI IMAGE -- around the peak
            I0roi = I0(ROIs{kroi}(1,1):ROIs{kroi}(2,1), ROIs{kroi}(1,2):ROIs{kroi}(2,2));
    
            %find Peak (max value)
            [maxval(kf),imax] = max(I0roi(:));
            [rmax,cmax] = ind2sub(size(I0roi),imax);

            scale=1.0;  %if(kf>2) scale = 1.5;   else scale=0.9; end;
            [level,EM] = graythresh(I0roi);
            I1 = im2bw(I0roi,level*scale);

            isel = find(I1(:)==1);
            SumCounts(kf,kroi) = sum(I0roi(isel));
            NumPix(kf,kroi) = length(isel);
            SumCountsBG(kf,kroi) = SumCounts(kf,kroi) - NumPix(kf,kroi)*BG_avg;


            %SHOW FIGURE
            [rsel,csel] = ind2sub(size(I0roi),isel);

            warning off

            if kroi==1          hf1 = figure;
                subplot(1,Nroi,1);
            else                figure(hf1);
                subplot(1,Nroi,kroi);
            end
            imshow(I0roi); hold all; plot(csel,rsel,'.','Color','r');
%         saveas(hf1,[toppath,'/',replace(filenames{kf},'.tiff','selectedPoints.jpg')], 'jpg');

        end
        %PRINT RESULT
        X = sprintf('%s , %d , %d, %d, %d, %d, %d, %d, %d',...
        filenames{kf},BG_avg, maxval_fullimg, ...
                        SumCounts(kf,1), SumCountsBG(kf,1), NumPix(kf,1),...
                        SumCounts(kf,2), SumCountsBG(kf,2), NumPix(kf,2) );
        disp(X);

        warning on;


    end





end


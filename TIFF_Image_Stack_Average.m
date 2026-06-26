% 1. Open stack of "identical" TIFF Images
% 2. Average them
 
%
function TIFF_Image_Stack_Avergae()
    close all;
    clear;
    program_name = 'image analysis';
    
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2024_11_04_LabData';
    start_path = 'G:\Shared drives\Corp Main\Engineering\LAB EXPERIMENTS\2025_02_10_LabData\P2P Round2 Measurements';
    
%% select IMAGE STACK
   
    %Select by Hand
    cd(start_path);
    [filenames,path] = uigetfile('*.PNG','Multiselect','on');
    Nf = length(filenames);
    if isstr(filenames)
        filename = filenames; clear filenames;
        filenames{1} = filename;
        Nf=1;
    end
 
%     [filename_1,path] = uigetfile('*.tiff','Multiselect','off');
%     [DAT.Iref] = get_image([path,filename_1]) ;
%     figure; imagesc(DAT.Iref)
%     figure; imagesc(Iavg);caxis([0,20])
    
    
%% Read All Files
    for kf = 1:Nf        
       [DAT.I0] = get_image([path,filenames{kf}]) ;                  
       if kf==1  %edefine sum-image
           Iavg = zeros(size(DAT.I0));           
       end
       Iavg = Iavg + DAT.I0;              
    end
    Iavg = Iavg / Nf;
    
    tmp = split(path,'\');
    title_name = replace(tmp{end-1},'_',' ');
    
    figure; imagesc(Iavg);  title(title_name);
    %caxis([0,20]);
    
    save([path,'/Iavg.mat'],'Iavg');
    imwrite(uint16(Iavg),[path,'\Iavg.tiff']);
    
    % Select ROI and cut
    hmsg = msgbox({'Select ROI. Two points: Upper Left and Lower Right'},'replace');
    set(hmsg,'Position',[300 300 250 150]);
    ROI=[0,0];
    [x0pick,y0pick]= ginput(1);  %grab pointer click from graph
    row1 = round(y0pick,0);  
    col1 = round(x0pick,0);
    ROI(1,:) = [row1,col1];
    hmsg2 = msgbox({'Upper Left selected.  Now Select Lower Right'},'replace');
    set(hmsg2,'Position',[300 300 250 150]);
    [x0pick,y0pick]= ginput(1);  %grab pointer click from graph
    row2 = round(y0pick,0);  
    col2 = round(x0pick,0);
    ROI(2,:) = [row2,col2];
  
 
    % compute histogram & average value & STDEV inside ROI
    [IMsel,meanVal,sumVal,Npix] = Analyze_ROI(ROI,Iavg,title_name);

    
    
    ROI2 = ROI + [100,0];
    
    % compute histogram & average value & STDEV inside ROI
    title_name = 'Background Level'
    [IMsel,meanVal,sumVal,Npix] = Analyze_ROI(ROI2,Iavg,title_name);
    
%     IMsel = Iavg(ROI(1,1):ROI(2,1),ROI(1,2):ROI(2,2));  %select ROI image
%     Vsel = IMsel(:);                                    %vector of pixel values
%     meanVal = mean(Vsel);
%     sumVal = sum(Vsel);
%     Npix = length(Vsel);
%     [yhist,xhist] = hist(Vsel,40);
%     figure; 
%     subplot(2,1,1); %imagesc(IMsel); colorbar;
%     imagesc(Iavg); colorbar;hold on; axis equal;
%     plot( [ROI(1,2),ROI(2,2),ROI(2,2),ROI(1,2),ROI(1,2)],...
%           [ROI(1,1),ROI(1,1),ROI(2,1),ROI(2,1),ROI(1,1)],'Color','k');
%     title({['MeanVal = ',num2str(meanVal)];...
%            ['SumVal = ',num2str(sumVal)];...
%            ['NPIX = ',num2str(Npix)];...
%             title_name});
%     subplot(2,1,2); bar(xhist,yhist); grid on;

    disp([num2str(meanVal),',',num2str(sumVal),',',num2str(Npix)]);
%     pause;

    
    
    
% %% select BG IMAGE
%     %Select by Hand
%     cd(start_path);
%     [BG_filename,BG_path] = uigetfile('*.tiff','Multiselect','off');
%     if BG_filename==0   BGflag = 0;
%     else                BFflag = 1;
%     end
%     
%     if isstr(BG_filename)   BG_filename = BG_filename;  
%     end
%     
%     
% %% Read BG Image
%     [BG.I0,BG.peak] = get_image([BG_path,BG_filename]) ;   
%  
% %% Define ROI
% ROI.r1 = BG.peak.r-50; 
% ROI.r2 = BG.peak.r+50; 
% ROI.c1 = BG.peak.c-50; 
% ROI.c2 = BG.peak.c+50; 

%% Generate Average
     

end  


%% Analyze ROI
function [IMsel,meanVal,sumVal,Npix] = Analyze_ROI(ROI,I0,title_name)
    % compute histogram & average value & STDEV inside ROI
    IMsel = I0(ROI(1,1):ROI(2,1),ROI(1,2):ROI(2,2));  %select ROI image
    Vsel = IMsel(:);                                    %vector of pixel values
    meanVal = mean(Vsel);
    sumVal = sum(Vsel);
    Npix = length(Vsel);
    [yhist,xhist] = hist(Vsel,40);
    figure; 
    subplot(2,1,1); %imagesc(IMsel); colorbar;
    imagesc(I0); colorbar;hold on; axis equal;
    plot( [ROI(1,2),ROI(2,2),ROI(2,2),ROI(1,2),ROI(1,2)],...
          [ROI(1,1),ROI(1,1),ROI(2,1),ROI(2,1),ROI(1,1)],'Color','k');
    title({['MeanVal = ',num2str(meanVal)];...
           ['SumVal = ',num2str(sumVal)];...
           ['NPIX = ',num2str(Npix)];...
            title_name});
    subplot(2,1,2); bar(xhist,yhist); grid on;

    disp([num2str(meanVal),',',num2str(sumVal),',',num2str(Npix)]);
end

%% Load TIFF Image + Find Center Peak

function [I0] = get_image(full_filename)
    I0  = imread(full_filename);
    if size(I0,3)==4
        I0  = I0 (:,:,1:3);
    end
    I0 = double(I0);

%     %convolve and find peak of iamge    
%     A = conv2(I0,ones(10))/(100);
%     [peak.maxval,imax] = max(A(:));
%     [peak.r,peak.c] = ind2sub(size(A),imax);
    
    
    s = split(full_filename,'\');
    filename = s(end);
    disp(['Got Image: ',filename]);
end


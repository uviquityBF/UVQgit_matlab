% Use MAT (NPZ)files from Overhead Loss to Generate HDR Images
 
function OverheadLoss_IntermediateData_HDR_Images_Sequence()
% Process overhead loss MAT files, detect defects, and export CSVs with metadata

    close all; clc;
    dirs.root = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run G\AM156\SD53\OverheadScattering_RunG_SD53_id4.5';

      
    HDRImages_Batch(dirs);

end

%% =========================================================
%  Main Batch Processing Function
%  =========================================================
function HDRImages_Batch(dirs)
    do_movie = 0;
    
    
    rootFolder = uigetdir(dirs.root, 'Select folder containing MAT files');
    if rootFolder == 0
        return;
    end
 
    %GET FILES in all siubfolders below...
    fileStruct = dir(fullfile(rootFolder, '**', '*.mat'));
    Nf = numel(fileStruct);

    for kf = 1:Nf

        path_and_file = fullfile(fileStruct(kf).folder, fileStruct(kf).name);
        %split path and file
        s = split(path_and_file,'\');
        filename = s{end};
        tmp = split(path_and_file,filename);
        path = tmp{1};

        %LOAD DATA
        load(path_and_file);
        
        % Package inputs
        hdrNormalized = double(hdrNormalized);
        a.input_facet   = input_facet;
        a.output_facet  = output_facet;
                
        %Z 
        s = split(path,'z='); s2 = split(s{end},'um'); 
        z_um(kf,1) = str2num(s2{1});

        %Image - generate movie
        if do_movie==1
            if kf==1 figure; end
            imagesc(log10(double(hdrNormalized))); colormap('bone');
            text(950,50,['z = ',num2str(z_um(kf,1)),' um'],'Color','w');
            title({['z=',num2str(z_um(kf,1))];replace(filename,'_',' ')});
            F(kf) = getframe;
        
        end
        
        %ROI
        ROIs(1).UL = [180,448];  %[x,y]  %horizontal, left
        ROIs(1).LR = [430,514];  %[x,y]
        label{1} = 'horizontal, left';
        ROIs(2).UL = [440,225];  %[x,y]  %vertical upper
        ROIs(2).LR = [537,493];  %[x,y]
        label{2} = 'vertical, upper';
        ROIs(3).UL = [555,404];  %[x,y]  %vertical center
        ROIs(3).LR = [660,612];  %[x,y]
        label{3} = 'vertical, center';
        ROIs(4).UL = [591,155];  %[x,y]  %horizontal, right
        ROIs(4).LR = [1175,225];  %[x,y]
        label{4} = 'horizontal, right';
        ROIs(5).UL = [923,655];  %[x,y]  %empty
        ROIs(5).LR = [1104,797];  %[x,y]
        label{5} = 'empty';
        
        %Get ROI Values
        for k=1:length(ROIs)
            Vals(kf,k) = processROI(hdrNormalized,ROIs(k));
        end
                
    end %loop over files
    
    
    if do_movie==1
        %playback
        movie(F,1,1);

        v = VideoWriter([dirs.root,'\newfile.avi'],'Uncompressed AVI');
        v.FrameRate = 2;
        open(v);
        writeVideo(v,F);
        close(v);
    end
    
    
    figure; 
    semilogy(z_um,Vals,'Linewidth',3); legend(label);
%     plot(z_um,Vals,'Linewidth',3); legend(label);
%     hp = plotyy(z_um,Vals(:,[1,2,4,5]),z_um,Vals(:,3)); legend(label);
    title('ROI signal levels v Z alignment position'); 
    
    figure; 
    semilogy(z_um,Vals,'.','Markersize',20); legend(label);

    title('ROI signal levels v Z alignment position'); 

end




function V = processROI(fullimage,ROI)
    I_ROI = fullimage(ROI.UL(2):ROI.LR(2)  ,  ROI.UL(1):ROI.LR(1) );
    
    %figure; imagesc(I_ROI);
  
    %Apply Threshold -->  Average only Values above threshold
    [yh,xh] = hist(I_ROI(:));
    t_frac = otsuthresh(yh);
    t = t_frac*(max(I_ROI(:))-min(I_ROI(:))) + min(I_ROI(:));
    isel = find(I_ROI(:)>t);
    V = mean(I_ROI(isel));    
    
    %Simple Avgerage of ALL PIXELS
    V = mean(I_ROI(:));
    
end

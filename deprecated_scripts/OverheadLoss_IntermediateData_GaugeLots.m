% Mess with NPZ files from Overhead Loss Output

% To do
% loop through all MAT files
% plot Intensity vector for each maximum deviation
%      - 
%
function OverheadLoss_IntermediateData_GaugeLots()
    clear; close all;

    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\GaugeLots\GL01\';
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run G\AM156\SD51';

    %Single File
    %OverheadLoss_IntermediateData_SingleFile(dirs);

    %Select Group of Files & Repeat over them   
    OverheadLoss_IntermediateData_Batch(dirs);
    
  
    
    
end

function OverheadLoss_IntermediateData_Batch(dirs)
    do_single_figure_plot = 1;

    rootFolder = uigetdir(dirs.d);
    fileStruct = dir(fullfile(rootFolder, '**', '*.mat'));      % GET LIST OF **ALL FILES**
    
    Nvmax = 0;    Nf = length(fileStruct);
    for kf = 1:Nf
        path_and_file = [fileStruct(kf).folder,'\',fileStruct(kf).name]; 
    
        %OverheadLoss_IntermediateData_SingleFile(path_and_file);

        %Open File
        matFileList{kf} = fullfile(fileStruct(kf).folder, fileStruct(kf).name);
        load([matFileList{kf}]); 

        %Get Defect List
        a.x_index_array = x_index_array;  a.y_data_array = y_data_array; 
        a.x_iqr = x_iqr;  a.y_iqr = y_iqr;  a.y_filtered = y_filtered; 
        a.fit_region_x = fit_region_x;  a.fit_region_y = fit_region_y;
        dd_list{kf} = Defect_Detection(a);                  %!!  GET DEFECTS  !!
        Nvmax = max([Nvmax,length(dd_list{kf}.istart)]);  %log how many values 
        
                
        %generate SINGLE FIGURE plot of all 
        if do_single_figure_plot==1
            Nplots = length(fileStruct);
            nr = floor(sqrt(Nplots));
            ncol = ceil(Nplots/nr);

            if kf==1 hf0=figure; else figure(hf0); end
            subplot(nr,ncol,kf); semilogy(x_index_array', y_data_array');
            hold all; semilogy(x_iqr,[y_iqr',y_filtered']);
            hold all; semilogy( fit_region_x,  fit_region_y,'.' );
            labelname = replace(replace(fileStruct(kf).name,'_',' '),'Missing','') ;
            labelname = replace(labelname,'intermediateData','');
            labelname = replace(labelname,'.mat','');
            xlabel( labelname );ylabel(['N defects = ',num2str(length(dd_list{kf}.x_um))]);
            file_labels{kf} = labelname;
        end
                       
    end
    
    %write defect data
    DDall.istart   = zeros(Nf,Nvmax);  
    DDall.iend     = zeros(Nf,Nvmax);  
    DDall.x_um     = zeros(Nf,Nvmax);  
    DDall.mag      = zeros(Nf,Nvmax);
    for kf=1:Nf
       Nd = length(dd_list{kf}.istart);
       DDall.istart(kf,1:Nd)    = dd_list{kf}.istart;
       DDall.iend(kf,1:Nd)      = dd_list{kf}.istart;
       DDall.x_um(kf,1:Nd)      = dd_list{kf}.x_um;
       DDall.mag(kf,1:Nd)       = dd_list{kf}.mag;       
    end
    disp(file_labels');
    csvwrite([dirs.d,'\DD_x_um.csv'],DDall.x_um);
    csvwrite([dirs.d,'\DD_mag.csv'],DDall.mag);
    
end



%% Defect Detection

function [dd] = Defect_Detection(a)
%     figure;  semilogy(a.x_index_array', a.y_data_array');
%     hold all; semilogy(a.x_iqr,[a.y_iqr',a.y_filtered']);
    
    y_iqr2 = interp1(a.x_iqr, double(a.y_iqr), a.x_index_array);
    D = ( a.y_data_array ~= y_iqr2);     %criterion:   y_iqr > y_data
    hw=5;

    D2 = merge1D(D,hw);  %merge together points within 'hw' pixels of each other
    dd.istart = find(diff(D2)==1);
    dd.iend = find(diff(D2)==-1);
    if dd.iend(1) < dd.istart(1)   dd.iend = dd.iend(2:end); end  %ensure that first value is not the 'end' of a defect domain
    %defect count
    Nd = min([length(dd.istart),length(dd.iend)]);
    dd.istart = dd.istart(1:Nd); dd.iend = dd.iend(1:Nd);
    dd.x_um = mean( [a.x_index_array(dd.istart); a.x_index_array(dd.iend)] ) ;
        
    %magnitude of defect
    delta_iqr = double(a.y_data_array) - y_iqr2;
    for k=1:Nd
        dd.mag(k,1) = sum(delta_iqr(dd.istart(k):dd.iend(k)));        
    end
    
	%debug plots
%     figure; plot(a.x_index_array,[D',D2']); ylim([-0.5,1.5]);
%     figure; stem(dd.x_um,ones(length(dd.x_um))); 
%     figure; plot(dd.x_um,ones(length(dd.x_um)),'.'); 
        
end

function x_merged = merge1D(x,hw)
    x_modify = zeros(size(x));
    N = length(x);
    for k=1:length(x)
        if ( x(k)==0 & sum(x(max([1,k-hw]):k)>0 & sum(x(k:min([k+hw],N)))>0 ) )
            x_modify(k) = 1;
        end
    end
    x_merged = x + x_modify;   
    
end


%% Get Data from all Files in Folder Structure
function remove_NaN()
        %remove NaN
        for kr=1:size(MAT.y_index,1)
            if sum(MAT.y_index)==NaN
            end
        end
end

function OverheadLoss_IntermediateData_SingleFile(dirs)

    %% Plot data from Single File
    if ~isfile(dirs)
        cd(dirs.d);
        [file,path] = uigetfile('*.mat');
        path_and_file = [path,file];
    else
        path_and_file = dirs;
        s = split(path_and_file,'\s');
        file = s{end};
    end
    load(path_and_file);
    
    
    figure;   semilogy(x_index_array', y_data_array');
    hold all; semilogy(x_iqr,[y_iqr',y_filtered']);
    hold all; semilogy( fit_region_x,  fit_region_y,'.' );
    label = replace(file,'Missing',''); label = replace(label,'_',' '); 
    title(label);
    
%     figure; imagesc(log10(double(wgSelectedImage)))
%     figure; imagesc(log10(double(rotatedCropped)))
    
    %% Run Fitting over many Domains
     % params = Fitting_Tests(x_iqr',y_iqr',1);
     
     
     %% Estimate BG 
      % BG = Estimate_BG(hdrNormalized,input_facet,output_facet,1);
%      theta_deg = atand( (output_facet(2)-input_facet(2)) / (output_facet(1)-input_facet(1)) );
%      hdrNormalized_horizontal = double(imrotate(hdrNormalized,theta_deg,'crop'));
%      figure; imagesc(log10(hdrNormalized_horizontal) );
%      c1 = uint16(input_facet(1))+100;
%      c2 = uint16(output_facet(1))-100;
%      V = sum(hdrNormalized_horizontal(:,c1:c2),2);
%      [dummy,imax] = max(V);
%      isel = [1:imax-80,imax+80:length(V)]';
%      BG = max(V(isel))
%      figure; plot(log10(V));
      
end

function BG = Estimate_BG(hdrNormalized,input_facet,output_facet,do_plot)
     theta_deg = atand( (output_facet(2)-input_facet(2)) / (output_facet(1)-input_facet(1)) );
     hdrNormalized_horizontal = double(imrotate(hdrNormalized,theta_deg,'crop'));
     c1 = uint16(input_facet(1))+100;
     c2 = uint16(output_facet(1))-100;
     V = mean(hdrNormalized_horizontal(:,c1:c2),2);
     [dummy,imax] = max(V);
     isel = [1:imax-80,imax+80:length(V)]';
     %BG = max(V(isel));
     Vsort = sort(V(isel),'descend');
     BG = mean(Vsort(1:50));
     if do_plot==1
        figure; imagesc(log10(hdrNormalized_horizontal) );
        figure; plot(log10(V));ylabel('log10(cts/pixel)'); grid on
     end
end

function params = Fitting_Tests(x,y,doplots)
    
    logY = log10(double(y));
    if doplots==1
        figure; plot(x,logY); hold on; ylabel('log10(y)'); xlabel('position, [um]');
    end
    
    xstart = [100:25:1000]';
    xend   = [1100:25:2000]';
    for ks = 1:length(xstart)
        for ke = 1:length(xend)
            isel = find( xstart(ks) <= x & x <= xend(ke) );
            tmp.p = polyfit(x(isel),logY(isel),1);
            slopes(ks,ke) = tmp.p(1);
            intercepts(ks,ke) = tmp.p(2);  
            if doplots==1
                plot(x(isel),polyval(tmp.p,x(isel)),'--','Color',[0.8,0.8,0.8]);
            end
        end   
    end
    Losses = 10*slopes * 1e4; % 10 * [log / um] * (10000um/cm)
    
    if doplots==1
        %Loss Values -- colormap
        figure; imagesc(xend, xstart,  Losses);
        xlabel('x end [um]');  ylabel('x start [um]'); title('Loss [dB/cm] calculated from many fit regions')
        colorbar;

        %Loss values -- histogram
        [yh,xh] = hist(Losses(:),20);
        figure; bar(xh,yh);
    end
    
    params = median(Losses(:));
    
    

end








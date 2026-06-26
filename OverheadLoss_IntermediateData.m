% Mess with NPZ files from Overhead Loss Output

% To do
% loop through all MAT files
% plot Intensity vector for each
% highlight & flag any potential failing criteria
% candidate criteria:
%      - maximum deviation
%      - 
%
function OverheadLoss_IntermediateData()
    clear; close all;

    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\GaugeLots\GL01\AM123\sd01\id3.7_vary_align_v_Z\Overhead Loss';
%    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\';

   % OverheadLoss_IntermediateData_SingleFile(dirs);
    
    OverheadLoss_IntermediateData_Batch_ManyFiles(dirs)
    
end

%% Get Data from all Files in Folder Structure
function OverheadLoss_IntermediateData_Batch_ManyFiles(dirs)
    do_corr = 1;
        
    %[0]
    rootFolder = dirs.d;
    fileStruct = dir(fullfile(rootFolder, '**', '*.mat'));  % GET LIST OF **ALL FILES**
    
    
    %[1] Identify group of files to compare -- based on attributes
    for k = 1:length(fileStruct)
        s = split(fileStruct(k).name,'_');
        tags.wfr{k}     =s{2};
        tags.subdie{k}  =s{3};
        tags.id{k}      =s{4};      
        sf = split(fileStruct(k).folder,'\');
        %tags.group{k}   = sf{end-1};
        s = split(sf{end},'_');
        tags.group{k}   = s{end};
    end
    tags.unique.id = unique(tags.id);
    tags.unique.subdie = unique(tags.subdie);
    tags.unique.wfr = unique(tags.wfr);
    tags.unique.group = unique(tags.group);
    
    %loop over each group (each device id)
    for ku = 1:length(tags.unique.id)
        this_tag = tags.unique.id{ku}
        k_sel = find(strcmp(tags.id,this_tag));
        group2id_list{ku} = tags.group{k_sel(1)};
        id_list{ku} = this_tag;
                
        %loop over members within a single group
        Nplots = length(k_sel);
        nr = floor(sqrt(Nplots));
        ncol = ceil(Nplots/nr);
        if Nplots<=36 do_plot=1;  else  do_plot=0; end
%         do_plot=1;
         
        for ik = 1:length(k_sel)  %loop over repetitions of same device
            k = k_sel(ik);
            source_this_instance{ik} = tags.group{k};
            matFileList{k} = fullfile(fileStruct(k).folder, fileStruct(k).name);
            load([matFileList{k}]);
            %generate plot of all 
            if do_plot==1
                if ik==1 hf0=figure; else figure(hf0); end
                subplot(nr,ncol,ik); semilogy(x_index_array', y_data_array');
                hold all; semilogy(x_iqr,[y_iqr',y_filtered']);
                hold all; semilogy( fit_region_x,  fit_region_y,'.' );
                xlabel([replace(fileStruct(k).name,'_',' ')]);
                ylabel(replace(source_this_instance{ik},'_',' '));
            end
            % generate data matrices on common X axes for correlatoin
            if ik==1
                commonX_index = double(x_index_array');
                commonX_iqr   = double(x_iqr');
                commonX_fit   = double(fit_region_x);                
            end
            MAT{ku}.y_index(:,ik)   = interp1(x_index_array',double(y_data_array'),commonX_index); 
            MAT{ku}.y_iqr(:,ik)     = interp1(x_iqr',double(y_iqr'),commonX_iqr); 
            MAT{ku}.y_filtered(:,ik)= interp1(x_iqr',double(y_filtered'),commonX_iqr); 
            MAT{ku}.y_fit(:,ik)     = interp1(fit_region_x,double(fit_region_y), commonX_fit);    
            
            %BG estimate
            BGs(ku,ik) = Estimate_BG(hdrNormalized,input_facet,output_facet,0);
            %Pout / Pin
            %inX = [100,300];
            inX = [900,1100];
            outX = [1900,2100];
            isel_in = find(inX(1) <= commonX_index & commonX_index <= inX(2) );
            isel_out = find(outX(1) <= commonX_index & commonX_index <= outX(2) );
            Throughput_Tot0(ku,ik) =  ( mean( MAT{ku}.y_index(isel_out,ik))) ...
                                   / ( mean( MAT{ku}.y_index(isel_in,ik))) ;
            Throughput_Tot(ku,ik) =  ( mean( MAT{ku}.y_index(isel_out,ik)) - BGs(ku,ik) ) ...
                                   / ( mean( MAT{ku}.y_index(isel_in,ik))- BGs(ku,ik) ) ;
            
            %SNR
            SNR =  (MAT{ku}.y_index(:,ik) - BGs(ku,ik)) ./ sqrt( MAT{ku}.y_index(:,ik) + BGs(ku,ik));             
            if do_plot==1  
                if ik==1 hf = figure; else figure(hf); end
                subplot(1,2,1); hold all; semilogy(commonX_index, MAT{ku}.y_index(:,ik) );  ylabel('value');
                susbplot(1,2,2); hold all; semilogy(commonX_index, SNR ); ylabel('SNR = (S-BG)/(S+BG)^0^.^5');
            end               
        end     
        
        %Correlations
        if do_corr==1
            CC{ku}.index        = corr(rmmissing(MAT{ku}.y_index));
            CC{ku}.iqr          = corr(rmmissing(MAT{ku}.y_iqr));
            CC{ku}.filtered     = corr(rmmissing(MAT{ku}.y_filtered));
            CC{ku}.fit          = corr(rmmissing(MAT{ku}.y_fit));        
            figure; 
            subplot(2,2,1); imagesc(CC{ku}.index);      xlabel('all data'); caxis([0.6,1]); colorbar; title(this_tag);
            subplot(2,2,2); imagesc(CC{ku}.iqr);        xlabel('iqr');      caxis([0.6,1]); colorbar;
            subplot(2,2,3); imagesc(CC{ku}.filtered);   xlabel('filtered'); caxis([0.6,1]); colorbar;
            subplot(2,2,4); imagesc(CC{ku}.fit);        xlabel('fit');      caxis([0.6,1]); colorbar;
        end
        
        %Average All Repeats of id='ku'
        if ku==1 commonX_all = commonX_index; end
        MAT_ALL1(:,ku) = interp1(commonX_index,mean(MAT{ku}.y_index,2),commonX_all); 
        MAT_ALL2(:,ku) = interp1(commonX_iqr,mean(MAT{ku}.y_iqr,2),commonX_all); 
        MAT_ALL3(:,ku) = interp1(commonX_iqr,mean(MAT{ku}.y_filtered,2),commonX_all); 
        MAT_ALL4(:,ku) = interp1(commonX_fit,mean(MAT{ku}.y_fit,2),commonX_all); 
        
        disp(source_this_instance');   
        figure; semilogy(commonX_index,mean(MAT{ku}.y_index,2));
        hold all; semilogy(commonX_iqr,[mean(MAT{ku}.y_iqr,2),mean(MAT{ku}.y_filtered,2)]);
        hold all; semilogy(commonX_fit,mean(MAT{ku}.y_fit,2),'.'); grid on;
        title([this_tag,' --  Average over all Repeats'])
        xlabel('position, x [um]');
        
    end  %repeat for each id    
    disp(group2id_list); size(group2id_list)
    
    disp(Throughput_Tot0)
    disp(BGs)
    disp(Throughput_Tot)
    
    figure; semilogy(commonX_all,MAT_ALL1);    legend(tags.unique.id);
    figure; semilogy(commonX_all,MAT_ALL2);    legend(tags.unique.id);
    figure; semilogy(commonX_all,MAT_ALL3);    legend(tags.unique.id);
    figure; semilogy(commonX_all,MAT_ALL4);    legend(tags.unique.id);
    
    % 
    id_exclude = {'id1.2';'id1.4';'id6.5';'id6.6'};
    ku_exclude = find(strcmp(id_exclude,id_exclude));    
    for kg = 1:length(tags.unique.group)
        ku_sel0 = find(strcmp(group2id_list, tags.unique.group{kg}));
        ku_sel = setdiff(ku_sel0,ku_exclude);
        disp(ku_sel);
        iref = find(250 <= commonX_all & commonX_all <= 650);
        YY1(:,kg) = mean(MAT_ALL1(:,ku_sel),2);     YY1(:,kg) = YY1(:,kg) / mean( YY1(iref,kg) );
        YY2(:,kg) = mean(MAT_ALL2(:,ku_sel),2);     YY2(:,kg) = YY2(:,kg) / mean( YY2(iref,kg) );
        YY3(:,kg) = mean(MAT_ALL3(:,ku_sel),2);     YY3(:,kg) = YY3(:,kg) / mean( YY3(iref,kg) );
        YY4(:,kg) = mean(MAT_ALL4(:,ku_sel),2);     YY4(:,kg) = YY4(:,kg) / max( YY4(:,kg) );            
    end
    figure; semilogy(commonX_all,YY1(:,1:4),'Linewidth',3); xlabel('X [um]');legend('250nm','350nm','700nm','1200nm'); %legend(tags.unique.group{1:4})
    figure; semilogy(commonX_all,YY2(:,1:4),'Linewidth',3); xlabel('X [um]');legend('250nm','350nm','700nm','1200nm'); %legend(tags.unique.group{1:4})
    figure; semilogy(commonX_all,YY3(:,1:4),'Linewidth',3); xlabel('X [um]');legend('250nm','350nm','700nm','1200nm'); %legend(tags.unique.group{1:4})
    figure; semilogy(commonX_all,YY4(:,1:4),'Linewidth',3); xlabel('X [um]');legend('250nm','350nm','700nm','1200nm'); %legend(tags.unique.group{1:4})
    
    figure; semilogy(commonX_all,YY1,'Linewidth',3); xlabel('X [um]');legend('250nm','350nm','700nm','1200nm','250-190','350,1deg','350,2deg','350,10deg'); %legend(tags.unique.group{1:4})
    figure; semilogy(commonX_all,YY2,'Linewidth',3); xlabel('X [um]');legend('250nm','350nm','700nm','1200nm','250-190','350,1deg','350,2deg','350,10deg'); %legend(tags.unique.group{1:4})
    figure; semilogy(commonX_all,YY3,'Linewidth',3); xlabel('X [um]');legend('250nm','350nm','700nm','1200nm','250-190','350,1deg','350,2deg','350,10deg'); %legend(tags.unique.group{1:4})
    figure; semilogy(commonX_all,YY4,'Linewidth',3); xlabel('X [um]');legend('250nm','350nm','700nm','1200nm','250-190','350,1deg','350,2deg','350,10deg'); %legend(tags.unique.group{1:4})
    
    
    
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [2] Multi-Domain Fitting on Every File 
    for k = 1:length(fileStruct)
        matFileList{k} = fullfile(fileStruct(k).folder, fileStruct(k).name);
        load([matFileList{k}]);
        %figure;  semilogy(x_index_array', y_data_array');
        %hold all; semilogy(x_iqr,[y_iqr',y_filtered']);
        %hold all; semilogy( fit_region_x,  fit_region_y,'.' );
        %title([replace(fileStruct(k).name,'_',' ')])
        
        params = Fitting_Tests(x_iqr',y_iqr',0);
        res{k,1} = tags.id{k};
        res{k,2} = tags.group{k};
        res{k,3} = params;
        
        disp(k);
    end
    disp(res{:,1});
    disp(res{:,2});
    disp(res{:,3});
    
    
end

function remove_NaN()
        %remove NaN
        for kr=1:size(MAT.y_index,1)
            if sum(MAT.y_index)==NaN
            end
        end
end

function OverheadLoss_IntermediateData_SingleFile(dirs)

    %% Plot data from Single File
    cd(dirs.d);
    [file,path] = uigetfile('*.mat');
    load([path,file]);
    
    figure;  semilogy(x_index_array', y_data_array');
    hold all; semilogy(x_iqr,[y_iqr',y_filtered']);
    hold all; semilogy( fit_region_x,  fit_region_y,'.' );
    
    figure; imagesc(log10(double(wgSelectedImage)))
    figure; imagesc(log10(double(rotatedCropped)))
    
    %% Run Fitting over many Domains
     params = Fitting_Tests(x_iqr',y_iqr',1);
     
     
     %% Estimate BG 
     BG = Estimate_BG(hdrNormalized,input_facet,output_facet,1);
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
%      
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








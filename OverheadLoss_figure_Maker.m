function OverheadLoss_figure_Maker()
    clear
    close all;
    
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\';
    filename = 'OverheadLoss_2025_12_08 - Sheet7.csv';
    
    %Read and Parse CSV data
    C = readCSV_headers_and_values([dirs.d,filename]);
    Run         = C(:,1);	
    wafer       = C(:,2);	
    SD          = C(:,3);
    id          = C(:,4);
    for kr=1:size(C,1)
        if length(str2num(C(kr,5)))>0   idn(kr,1)   = str2num(C(kr,5)); end
        if length(str2num(C(kr,6)))>0   Loss(kr,1)  = str2num(C(kr,6)); end
        if length(str2num(C(kr,7)))>0   Wd(kr,1)    = str2num(C(kr,7)); end
        if length(str2num(C(kr,8)))>0   Wa(kr,1)    = str2num(C(kr,8)); end
    end     
    
%     figure; plot(Wd,Loss,'.')
    
    
    %[1]  Remove Loss < 10dB/cm  
    minLoss_to_keep = 15;    
    isel        = find(Loss > minLoss_to_keep | Wd > 500);
    Run         = Run(isel);
    wafer       = wafer(isel);
    SD          = SD(isel);
    id          = id(isel);
    idn         = idn(isel);
    Loss        = Loss(isel);
    Wd          = Wd(isel);
    Wa          = Wa(isel);
    
    
    %[2] Group by Run,Wfr,SD
    isel_g{1}   = find(strcmp(Run,"RunA"));     grp_name{1}='RunA';
    isel_g{2}   = find(strcmp(Run,"RunB"));     grp_name{2}='RunB';
    isel_g{3}   = find(strcmp(Run,"Round2"));   grp_name{3}='Run2';
    isel_g{4}   = find(strcmp(wafer,"AM123"));  grp_name{4}='RunC 123';
    isel_g{5}   = find(strcmp(wafer,"KW"));     grp_name{5}='RunC KW';
    isel_g{6}   = find(strcmp(wafer,"AM135"));  grp_name{6}='RunF 135';
    isel_g{7}   = find(strcmp(wafer,"AM136"));  grp_name{7}='RunF 136';
    n=0; for k=1:7  n = n+length(isel_g{k}); end  
    disp(n);
    
    
    
    %[3]  For Each Group...
    Ng = length(isel_g);
    hf = figure;
    for kg=1:Ng
        disp(grp_name{kg});
        
        %select members of the group
        isel = isel_g{kg};
        tmp.Loss    = Loss(isel);
        tmp.Wd      = Wd(isel);
        tmp.Wa      = Wa(isel);
        tmp.idn     = idn(isel);
        
        %eliminate repeats for each device ID
        uID{kg}         = unique(tmp.idn);
        Nu              = length(uID{kg});
            %loop over every unique ID and boil down to single Loss value (Median)
            for ku=1:Nu  
                iselu = find(tmp.idn==uID{kg}(ku));
                uLoss{kg}(ku,1) = median( tmp.Loss(iselu) );
                uWd{kg}(ku,1) = median( tmp.Wd(iselu) );
                uWa{kg}(ku,1) = median( tmp.Wa(iselu) );
            end        
        
        %Average over Width Bin
        width_bin_full = 50;     %nm
        if (length(uWd{kg})>0)
            Wd_bins = [min(Wd):width_bin_full:max(Wd)+width_bin_full]';
            for kb = 1:length(Wd_bins)
                isel = find(  Wd_bins(kb)-0.5*width_bin_full  <= uWd{kg} ...
                            &  uWd{kg}<= Wd_bins(kb)+0.5*width_bin_full );
                Loss_v_Wd(kb,kg) = median(uLoss{kg}(isel));   
                LossSD_v_Wd(kb,kg) = std(uLoss{kg}(isel));                
                LossRange_v_Wd(kb,kg) = max([NaN,max(uLoss{kg}(isel))-min(uLoss{kg}(isel))]);                
            end
            Wa_bins = [min(Wa):width_bin_full:max(Wa)+width_bin_full]';    
            for kb = 1:length(Wa_bins)
                isel = find(  Wa_bins(kb)-0.5*width_bin_full  <= uWa{kg} ...
                            & uWa{kg} <= Wa_bins(kb)+0.5*width_bin_full );
                Loss_v_Wa(kb,kg) = median(uLoss{kg}(isel));      
                LossSD_v_Wa(kb,kg) = std(uLoss{kg}(isel));      
                LossRange_v_Wa(kb,kg) = max([NaN,max(uLoss{kg}(isel))-min(uLoss{kg}(isel))]);                
            end
        end
        
        %add data set to scattter plot
        figure(hf);
        plot(uWd{kg},uLoss{kg},'.','MarkerSize',15); hold all;
  
        
    end %loop over groups
    legend(grp_name); title('All Devices');

    

 %%  Figures
    hf2 = figure;  plot(Wd_bins,Loss_v_Wd,'.','MarkerSize',20); xlabel('width [nm] (drawn)'); grid on;   legend(grp_name);
    hf3 = figure;  plot(Wa_bins,Loss_v_Wa,'.','MarkerSize',20); xlabel('width [nm] (actual)'); grid on;  legend(grp_name);  
    
    
    hf4 = figure;  
    for kg=1:Ng
        scatter(Wa_bins, Loss_v_Wa(:,kg), 'filled'); grid on;   hold on;
        eb(1) = errorbar(Wa_bins, Loss_v_Wa(:,kg), LossSD_v_Wa(:,kg), 'vertical', 'LineStyle', 'none');
        %eb(2) = errorbar(Wa_bins, Loss_v_Wa, 0, 'horizontal', 'LineStyle', 'none');
        set(eb);
    end
    ylim([0,100]); legend(grp_name);  
    
    %% write output
    out1 = [Wd_bins,Loss_v_Wd,LossSD_v_Wd,LossRange_v_Wd];
    out2 = [Wa_bins,Loss_v_Wa,LossSD_v_Wa,LossRange_v_Wa];
    csvwrite([dirs.d,'\tmp.csv'],[out1;out2]);
    
     
end


%% ========================================================================
function C = readCSV_headers_and_values(filepath)
    % Read CSV using detectImportOptions, convert to strings
    opts = detectImportOptions(filepath);
    T = readtable(filepath, opts);
    C = string(table2cell(T));

    headers = C(1,:);    % first row
    values  = C(2:end,:);    % rest
end

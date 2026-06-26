% Use MAT (NPZ)files from Overhead Loss to DETECT DEFECTS

% Refactored with CHATGPT

function OverheadLoss_IntermediateData_DEFECT_ANALYSIS()
% Process overhead loss MAT files, detect defects, and export CSVs with metadata

    close all; clc;

%  dirs.root = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\GaugeLots\GL00\am123\sd01\Standardized Testing';
%  dirs.root = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\GaugeLots\GL01\AM125\sd01\StandardizedTesting - copyBF'
%  dirs.root = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\GaugeLots\GL01\AM125\sd01\Standardized Testing';
%  dirs.root = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\GaugeLots\GL01\AM125\sd01\StandardizedTesting - copyBF\Overhead Loss\Post Clean\';
%  dirs.root = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\Run G\AM156\SD51';
%  dirs.root = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\GaugeLots\GL02\AM123\Standardized Testing\Overhead Loss';
dirs.root = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\GaugeLots\GL01\AM125\sd01\Standardized Testing\Overhead Loss\Post Clean.2\'
    doSingleFigurePlot  = false;
    mergeHalfWidth      = 10;   % pixels for merging adjacent defects
    facet_exclusion_um  = 300;  %[um] - how far away the defects must be from the facets to be counted
    thresh              = 3;  
    
    OverheadLoss_DefectDetection_Batch(dirs,mergeHalfWidth,thresh,facet_exclusion_um,doSingleFigurePlot);

end

%% =========================================================
%  Main Batch Processing Function
%  =========================================================
function OverheadLoss_DefectDetection_Batch(dirs,mergeHalfWidth,thresh,facet_exclusion_um,doSingleFigurePlot)

    rootFolder = uigetdir(dirs.root, 'Select folder containing MAT files');
    if rootFolder == 0
        return;
    end
    
    % DELETE EXISTING PNG FILES ?
%     fileStruct_PNG = dir(fullfile(rootFolder, '**', '*DefectIDs*.png'))
%     for k=1:length(fileStruct_PNG)
%         delete([fileStruct_PNG(k).folder,'\',fileStruct_PNG(k).name]);
%         disp(['Deleted:  ',fileStruct_PNG(k).folder,'\',fileStruct_PNG(k).name]);
%     end
    
    fileStruct = dir(fullfile(rootFolder, '**', '*.mat'));
    Nf = numel(fileStruct);

    dd_list = cell(Nf,1);
    file_labels = cell(Nf,1);
    Nvmax = 0;

    % Metadata storage
    meta(Nf) = struct( ...
        'Run',"", 'Lot',"", 'Sample',"", ...
        'Waveguide',"", 'Date',"");

    if doSingleFigurePlot
        hf = figure;
    end

    % =====================
    % Loop over files
    % =====================
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
        a.x_index_array = x_index_array;
        a.y_data_array  = y_data_array;
        a.x_iqr         = x_iqr;
        a.y_iqr         = y_iqr;
        a.y_filtered    = y_filtered;
        a.fit_region_x  = fit_region_x;
        a.fit_region_y  = fit_region_y;
        a.input_facet   = input_facet;
        a.output_facet  = output_facet;
        

        % Metadata
        meta(kf) = parseAllIdentifiers(path_and_file);

        % Defect detection
        dd_list{kf} = Defect_Detection(a, mergeHalfWidth, thresh, facet_exclusion_um);
        Nvmax = max(Nvmax, numel(dd_list{kf}.istart));
        
        % Plot
        if doSingleFigurePlot
            subplotLayout(Nf, kf);
            semilogy(a.x_index_array, a.y_data_array, 'Color',[0.6 0.6 0.6]); hold on;
            % IQR curve (downsampled)
            semilogy(a.x_iqr, a.y_iqr, 'LineWidth', 1.2);
            % Filtered data (full resolution)
            semilogy(a.x_iqr, a.y_filtered, 'LineWidth', 1.2);
            % Fit region
            semilogy(a.fit_region_x, a.fit_region_y, '.', 'MarkerSize', 6);
            hold off;

            label = cleanFileLabel(fileStruct(kf).name);
            xlabel(label);
            ylabel(sprintf('N defects = %d', numel(dd_list{kf}.x_um)));
            file_labels{kf} = label;
        end
        
        %Show Defect IDs on Image
        dd_list{kf}.x_pix = mean([dd_list{kf}.istart;dd_list{kf}.iend],1) + input_facet(1);
        dd_list{kf}.y_pix = mean([input_facet(2),output_facet(2)])*ones(size(dd_list{kf}.x_um));
        y1=input_facet(2); y2=output_facet(2);
        dd_list{kf}.y_pix = y1 + (y2-y1)/(output_facet(1)-input_facet(1)) * (dd_list{kf}.x_pix - input_facet(1));
        
        sz = 8*abs(dd_list{kf}.mag);
        iselnan = find(isnan(sz)|sz==0); sz(iselnan)=0.1;
        
        hf=figure; imagesc(log10(double(hdrNormalized))); colormap(bone);
        hold on; scatter(dd_list{kf}.x_pix,dd_list{kf}.y_pix,sz,'r')
        ylim([input_facet(2)-100,input_facet(2)+100]); axis equal;
        title({'Defect IDs: ';replace(filename,'_',' ')});
        path=tmp{1};
        saveas(gcf,[path,filename(1:end-4),'_DefectIDs_thresh=',num2str(thresh),'.png']);
        close(hf);
%         hold on; plot(dd_list{kf}.x_pix,dd_list{kf}.y_pix+15,'^','Color','r');        
    end

    %% =====================================================
    %  Wide-format outputs (legacy compatible)
    %  =====================================================
    DDall.istart = zeros(Nf, Nvmax);
    DDall.iend   = zeros(Nf, Nvmax);
    DDall.x_um   = zeros(Nf, Nvmax);
    DDall.mag    = zeros(Nf, Nvmax);

    for kf = 1:Nf
        Nd = numel(dd_list{kf}.istart);
        DDall.istart(kf,1:Nd) = dd_list{kf}.istart;
        DDall.iend(kf,1:Nd)   = dd_list{kf}.iend;
        DDall.x_um(kf,1:Nd)   = dd_list{kf}.x_um;
        DDall.mag(kf,1:Nd)    = dd_list{kf}.mag;
    end

    T_x = array2table(DDall.x_um);
    T_mag = array2table(DDall.mag);

    T_x = addvars(T_x, {meta.Run}', {meta.Lot}', {meta.Sample}', ...
                  {meta.Waveguide}', {meta.Date}', ...
                  'Before',1, ...
                  'NewVariableNames',{'Run','Lot','Sample','Waveguide','Date'});

    T_mag = addvars(T_mag, {meta.Run}', {meta.Lot}', {meta.Sample}', ...
                  {meta.Waveguide}', {meta.Date}', ...
                  'Before',1, ...
                  'NewVariableNames',{'Run','Lot','Sample','Waveguide','Date'});

    writetable(T_x,  fullfile(dirs.root,'DD_x_um.csv'));
    writetable(T_mag,fullfile(dirs.root,'DD_mag.csv'));

    %% =====================================================
    %  Master defect-level table
    %  =====================================================
    Master = table();
    Summary = table();
    for kf = 1:Nf
        dd = dd_list{kf};
        Nd = numel(dd.x_um);
        if Nd == 0, continue; end

        T = table();
        % %         T.Run         = repmat([pad(meta(kf).Run,100)], Nd, 1);
        %         T.Lot         = repmat([meta(kf).Lot,' '], Nd, 1);
        %         T.Sample      = repmat([meta(kf).Sample,' '], Nd, 1);
        %         T.Waveguide   = repmat([meta(kf).Waveguide,' '], Nd, 1);
        %         T.Date        = repmat([meta(kf).Date,' '], Nd, 1);
        T.Run       = repmat(string(meta(kf).Run), Nd, 1);
        T.Lot       = repmat(string(meta(kf).Lot), Nd, 1);
        T.Sample    = repmat(string(meta(kf).Sample), Nd, 1);
        T.Waveguide = repmat(string(meta(kf).Waveguide), Nd, 1);
        T.Date      = repmat(string(meta(kf).Date), Nd, 1);
        T.DefectIndex = (1:Nd).';
        T.x_um        = dd.x_um(:);
        T.mag         = dd.mag(:);
                
        Master = [Master; T]; %#ok<AGROW>
        
        %Summary Stats for this WG ID
        S = table();
        S.Run           = string(meta(kf).Run);
        S.Lot           = string(meta(kf).Lot);
        S.Sample        = string(meta(kf).Sample);
        S.Waveguide     = string(meta(kf).Waveguide);
        S.Nd            = Nd;
        S.MagSum        = sum(dd.mag);
        S.MagAvg        = mean(dd.mag);
        Summary = [Summary;S];
    end

     writetable(Summary,fullfile(dirs.root,'DD_Summary.csv'));
     writetable(Master, fullfile(dirs.root,'DD_Master_AllDefects.csv'));
    
    %% =====================================================
    %  Scatter plot of defect locations vs waveguide
    %  =====================================================

    % Extract numeric waveguide index
    WGnum = extractWaveguideNumber(Master.Waveguide);

    % Remove rows where WG number could not be parsed
    valid = ~isnan(WGnum);
    X = Master.x_um(valid);
    Y = WGnum(valid);
    M = Master.mag(valid);

    %
    isel = find(X>0);
    Msel_sort = sort(M(isel),'descend');
    M_max = median(Msel_sort(1:5))  %scale to "near-max" value
    
    % Sort by waveguide number
    [Ysorted, order] = sort(Y);
    Xsorted = X(order(isel));
    Msorted = M(order(isel));
    Msorted_size = sqrt(Msorted)/sqrt(M_max)*40;

     
    figure('Position',[500,50,575,1100]);
    subplot(6,1,1:5)
    scatter(Xsorted, Ysorted(isel), ...
            Msorted_size+0.1, log10(Msorted), 'filled');   % size=25, color=mag
    colormap(flipud(bone));
    c = colorbar;
    c.Label.String = 'Log10(Defect Magnitude)';
    caxis([0,2]);

    xlabel('X Position (\mum)');
    ylabel('Waveguide ID (numeric)');
    s = split(rootFolder,'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples');
    title({replace(s{end},'_',' ');'Defect Map Across Waveguides'});

    grid on;
    set(gca,'YDir','normal');

    subplot(6,1,6)
    [yh,xh] = hist(log10(Msorted),20); bar(xh,yh); xlabel('log10(magnitude)'); ylabel('number of defects');
    text(mean([median(xh),min(xh)]),max(yh)/2, 'Defect Count for all data' );
    
    saveas(gcf,[rootFolder,'\DefectMap_thresh=',num2str(thresh),'.png']);
    
  
    
end

%% =========================================================
%  Defect Detection
%  =========================================================
function dd = Defect_Detection(a, hw, thresh, facet_exclusion_um)

    y_iqr_interp = interp1(a.x_iqr, double(a.y_iqr), ...
        a.x_index_array, 'linear', 'extrap');
    
    defectMask = ( abs( double(a.y_data_array) - y_iqr_interp) ./ y_iqr_interp  >= thresh ) ;
    defectMask = merge1D(defectMask, hw);

    %debugging
%     figure; plot(a.x_index_array',[double(a.y_data_array'), y_iqr_interp'])
%     testvalues = abs( double(a.y_data_array) - y_iqr_interp) ./ y_iqr_interp;
%     figure; plot(testvalues);
    
    
    dd.istart = find(diff([0 defectMask]) == 1);
    dd.iend   = find(diff([defectMask 0]) == -1);

    Nd = min(numel(dd.istart), numel(dd.iend));
    dd.istart = dd.istart(1:Nd);
    dd.iend   = dd.iend(1:Nd);

    dd.x_um = mean([a.x_index_array(dd.istart); ...
                    a.x_index_array(dd.iend)], 1);

    %REMOVE DEFECTS NEAR FACETS
    ikeep = find( abs(dd.x_um - a.x_index_array(1)) > facet_exclusion_um  &...
                    abs(dd.x_um - a.x_index_array(end)) > facet_exclusion_um  )
    dd.x_um = dd.x_um(ikeep);
    dd.istart = dd.istart(ikeep);
    dd.iend = dd.iend(ikeep);
    Nd = length(ikeep);
    
                
    % DEFINITION OF MAGNITUDE
    delta = double(a.y_data_array) - y_iqr_interp;
    dd.mag = zeros(Nd,1);
    for k = 1:Nd
        dd.mag(k) = abs(  sum(delta(dd.istart(k):dd.iend(k)))  /  mean(y_iqr_interp(dd.istart(k):dd.iend(k)))  ) ;
    end

end

%% =========================================================
%  Merge Adjacent Defects
%  =========================================================
function xMerged = merge1D(x, hw)

    xMerged = x;
    N = numel(x);

    for k = 1:N
        if x(k) == 0
            left  = max(1, k-hw);
            right = min(N, k+hw);
            if any(x(left:k-1)) && any(x(k+1:right))
                xMerged(k) = 1;
            end
        end
    end
end

%% =========================================================
%  Metadata Parsing
%  =========================================================
function meta = parseAllIdentifiers(path_and_file)

    meta.Run       = "UnknownRun";
    meta.Lot       = "UnknownLot";
    meta.Sample    = "UnknownSample";
    meta.Waveguide = "UnknownWG";
    meta.Date      = "";

    parts = split(path_and_file, filesep);
    fname = parts{end};

    for k = 1:numel(parts)
        token = parts{k};

        if regexp(token,'^(Run|GL)\w*','once')
            meta.Run = token;
        end
        if regexp(token,'^AM\d+','once')
            meta.Lot = token;
        end
        if regexp(token,'^SD\d+','once')
            meta.Sample = token;
        end
        if regexp(token,'^(WG\d+|id\d+(\.\d+)?)','once')
            meta.Waveguide = token;
        end
    end

    if meta.Waveguide == "UnknownWG"
        m = regexp(fname,'(WG\d+|id\d+(\.\d+)?)','match','once');
        if ~isempty(m), meta.Waveguide = m; end
    end

    d = regexp(path_and_file,'\d{4}[-_]\d{2}[-_]\d{2}','match','once');
    if ~isempty(d)
        meta.Date = strrep(d,'_','-');
    end

    meta.Waveguide = normalizeWaveguide(meta.Waveguide);

end

function wg = normalizeWaveguide(wg)
    wg = string(wg);
    if startsWith(wg,"WG")
        n = regexp(wg,'\d+','match','once');
        wg = "WG" + string(str2double(n));
    end
end

%% =========================================================
%  Plot helpers
%  =========================================================
function subplotLayout(N, idx)
    nr = floor(sqrt(N));
    nc = ceil(N / nr);
    subplot(nr, nc, idx);
end

function label = cleanFileLabel(fname)
    label = erase(fname, {'_', 'Missing', 'intermediateData', '.mat'});
end

function wgnum = extractWaveguideNumber(wg)
% Convert WG identifiers to numeric values for plotting
% Examples:
%   WG1    -> 1
%   WG01   -> 1
%   id4.1  -> 4.1

    wg = string(wg);
    wgnum = nan(size(wg));

    for k = 1:numel(wg)
        token = regexp(wg(k), '\d+(\.\d+)?', 'match', 'once');
        if ~isempty(token)
            wgnum(k) = str2double(token);
        end
    end
end

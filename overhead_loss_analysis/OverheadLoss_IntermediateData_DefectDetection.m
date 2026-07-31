% Use MAT (NPZ)files from Overhead Loss to DETECT DEFECTS

% Refactored with CHATGPT

function OverheadLoss_IntermediateData_DEFECT_ANALYSIS()
% Process overhead loss MAT files, detect defects, and export CSVs with metadata

    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\overhead_loss');
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
    defect_criterion    = 'percent_deviation';  % 'percent_deviation' (default) or 'gauge_lots' (hw=5 exact-inequality, from the merged GaugeLots script)

    OverheadLoss_DefectDetection_Batch(dirs,mergeHalfWidth,thresh,facet_exclusion_um,doSingleFigurePlot,defect_criterion);

end

%% =========================================================
%  Main Batch Processing Function
%  =========================================================
function OverheadLoss_DefectDetection_Batch(dirs,mergeHalfWidth,thresh,facet_exclusion_um,doSingleFigurePlot,defect_criterion)

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
    
    [fileStruct, filePaths, matData] = load_intermediate_mats(rootFolder);
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

        path_and_file = filePaths{kf};
        filename = fileStruct(kf).name;
        path = [fileStruct(kf).folder, filesep];

        % Package inputs
        D = matData{kf};
        a.x_index_array = D.x_index_array;
        a.y_data_array  = D.y_data_array;
        a.x_iqr         = D.x_iqr;
        a.y_iqr         = D.y_iqr;
        a.y_filtered    = D.y_filtered;
        a.fit_region_x  = D.fit_region_x;
        a.fit_region_y  = D.fit_region_y;
        a.input_facet   = D.input_facet;
        a.output_facet  = D.output_facet;

        % Metadata
        meta(kf) = parse_lab_identifiers(path_and_file);

        % Defect detection
        dd_list{kf} = Defect_Detection(a, mergeHalfWidth, thresh, facet_exclusion_um, defect_criterion);
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
        dd_list{kf}.x_pix = mean([dd_list{kf}.istart;dd_list{kf}.iend],1) + a.input_facet(1);
        dd_list{kf}.y_pix = mean([a.input_facet(2),a.output_facet(2)])*ones(size(dd_list{kf}.x_um));
        y1=a.input_facet(2); y2=a.output_facet(2);
        dd_list{kf}.y_pix = y1 + (y2-y1)/(a.output_facet(1)-a.input_facet(1)) * (dd_list{kf}.x_pix - a.input_facet(1));

        sz = 8*abs(dd_list{kf}.mag);
        iselnan = find(isnan(sz)|sz==0); sz(iselnan)=0.1;

        hf=figure; imagesc(log10(double(D.hdrNormalized))); colormap(bone);
        hold on; scatter(dd_list{kf}.x_pix,dd_list{kf}.y_pix,sz,'r')
        ylim([a.input_facet(2)-100,a.input_facet(2)+100]); axis equal;
        title({'Defect IDs: ';replace(filename,'_',' ')});
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
function dd = Defect_Detection(a, hw, thresh, facet_exclusion_um, defect_criterion)
% defect_criterion selects how a per-pixel defect candidate is flagged:
%   'percent_deviation' (default) - |signal - iqr_fit| / iqr_fit >= thresh
%   'gauge_lots'                  - signal ~= iqr_fit (exact-inequality test,
%                                   merged in from OverheadLoss_IntermediateData_GaugeLots.m,
%                                   which always used hw=5)
% Everything downstream (merge width, facet exclusion, magnitude calc) is
% shared between both modes.

    if nargin < 5 || isempty(defect_criterion)
        defect_criterion = 'percent_deviation';
    end

    y_iqr_interp = interp1(a.x_iqr, double(a.y_iqr), ...
        a.x_index_array, 'linear', 'extrap');

    switch defect_criterion
        case 'percent_deviation'
            defectMask = ( abs( double(a.y_data_array) - y_iqr_interp) ./ y_iqr_interp  >= thresh ) ;
        case 'gauge_lots'
            defectMask = ( double(a.y_data_array) ~= y_iqr_interp );
        otherwise
            error('Defect_Detection:UnknownCriterion','Unknown defect_criterion "%s"',defect_criterion);
    end
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

% Metadata parsing (Run/Lot/Sample/Waveguide/Date) now lives in the shared
% helpers\overhead_loss\parse_lab_identifiers.m -- see call in
% OverheadLoss_DefectDetection_Batch above.

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

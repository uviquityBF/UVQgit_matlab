% Generate Surface for index (of a mode) = f( X1, X2) 
% Input format assumed (for surface definition):      
%     - first  Col: X1 values (WLs)
%     - second Col: X2 values (widths)
%     - third Col:  Y values  (index)
%
%

function surface_comparison_forIndexMatching()
    clear; close all
    
    start_path = 'G:\Shared drives\Corp Main\Engineering\Design & SHG Analysis (Brent)\COMSOL Modeling\Results\modal dispersion tables (for matlab)';
    WLs     = [200:5:400];      %WL  (SHG)
    widths  = [140:20:2000];    %width    
    WLs     = [210:2:300];      %WL  (SHG)
    widths  = [140:20:800];     %width    
    
    dN_f = +0;
    dN_s = +0;

    %% GET SURFACES 
    [F.WL,  F.widths,   F.vq, F.infotext]          = get_surface_definition(start_path, WLs, widths, dN_f);   %% n_FUND = f(WL,width)

    %[SHG20.WL, SHG.widths, SHG.vq, SHG.infotext]     = get_surface_definition(start_path, WLs, widths, dN_s);   %% n_shg = f(WL,width)
    [SHG20.WL, SHG20.widths, SHG20.vq, SHG20.infotext]     = get_surface_definition(start_path, WLs, widths, dN_s);   %% n_shg = f(WL,width)
    [SHG40.WL, SHG40.widths, SHG40.vq, SHG40.infotext]     = get_surface_definition(start_path, WLs, widths, dN_s);   %% n_shg = f(WL,width)
    [SHG04.WL, SHG04.widths, SHG04.vq, SHG04.infotext]     = get_surface_definition(start_path, WLs, widths, dN_s);   %% n_shg = f(WL,width)

    
    width2plt = 300;
    v00 = interp2(F.WL, F.widths, F.vq, WLs,width2plt*ones(size(WLs)) )';
    v20 = interp2(SHG20.WL, SHG20.widths, SHG20.vq, WLs,width2plt*ones(size(WLs)) )';
    v40 = interp2(SHG40.WL, SHG40.widths, SHG40.vq, WLs,width2plt*ones(size(WLs)) )';
    v04 = interp2(SHG04.WL, SHG04.widths, SHG04.vq, WLs,width2plt*ones(size(WLs)) )';
    figure; plot(WLs,[v00,v20,v40,v04],'Linewidth',2); ylim([1.9,2.2]); xlabel('Wavelength [nm] (of SH)');  ylabel('Effective Index of Mode');  
    legend('TM_0_0(f)','TM_2_0','TM_4_0','TM_0_4');  grid on; title(['width = ',num2str(width2plt),'nm']);
    
    
	WL2plt = 226;
    v00w = interp2(F.WL, F.widths, F.vq,  WL2plt*ones(size(widths)), widths )';
    v20w = interp2(SHG20.WL, SHG20.widths, SHG20.vq, WL2plt*ones(size(widths)), widths)';
    v40w = interp2(SHG40.WL, SHG40.widths, SHG40.vq,  WL2plt*ones(size(widths)), widths )';
    v04w = interp2(SHG04.WL, SHG04.widths, SHG04.vq,  WL2plt*ones(size(widths)), widths )';
    figure; plot(widths,[v00w,v20w,v40w,v04w],'Linewidth',2); ylim([1.9,2.2]);  xlabel('Waveguide Width [nm]');  ylabel('Effective Index of Mode');  
    legend('TM_0_0(f)','TM_2_0','TM_4_0','TM_0_4');  grid on; title(['WL = ',num2str(WL2plt),'nm']);
    xlim([150,800]);

    
    titletext = replace([F.infotext, ' + ', SHG.infotext],'temp surfacedefinition',' ');
    
    %plot both surfaces over each other
    hf(1) = figure; 
    mesh(F.WL,     F.widths,   F.vq);    hold on;
    mesh(SHG.WL, SHG.widths, SHG.vq); 
    xlabel('WL');ylabel('width');  
    title(titletext);
    
    
    %% NUMERICAL:   MESHGRID + Numerical Intersection
    %Find Surface Intersections - NUMERICAL
	[tols,minskeep,fine] = find_surface_intersections_numerical(F,SHG,titletext);
    
    
    %% ANALYTICAL:  PolyFits + Analytic Surfaces
    
    %reference point for Taylor Expansions
    X0.WL       = 230;
    X0.width    = 350;
    
    %polyfit v. WL  **at each WIDTH** (row)
    [F.pp_v_WL]      = polyfit_manycurves_byrow(F.WL(1,:) ,   F.vq,   X0.WL, 'FUND:  fit v. WAVELENGTH[nm]');
    [SHG.pp_v_WL]    = polyfit_manycurves_byrow(SHG.WL(1,:) , SHG.vq, X0.WL, 'SHG:  fit v. WAVELENGTH[nm]');
    
    %polyfit v. WIDTH  **at each WL** (col)
    infotext='fit v. width[nm]';
    [F.pp_v_wdt,]      = polyfit_manycurves_byrow(F.widths(:,1)'   , F.vq',    X0.width, 'FUND:  fit v. WIDTH[nm]')';
    [SHG.pp_v_wdt]    = polyfit_manycurves_byrow(SHG.widths(:,1)' , SHG.vq',  X0.width,  'SHG:  fit v. WIDTH[nm]')';
    
    %plot Fit-surfaces
    vec_WL          = WLs';     % [200:10:500]';
    vec_wdt         = widths';  % [300:10:700]';
    [WL_roots] = find_surface_intersections_analytic(F, SHG, X0, vec_WL, vec_wdt);
    titletext = {[replace(F.infotext,'surfacedefinition',' ')  ,'+',  replace(SHG.infotext,'surfacedefinition',' ') ];'Analytic Intersection'};
    title(titletext);
    
    disp([WL_roots.pos, WL_roots.neg, vec_wdt]);
     
    %% Save all Open Figures
    FolderName = uigetdir(start_path)
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      %FigName   = get(FigHandle, 'Name');
      FigName = [num2str(iFig),'.fig'];
      savefig(FigHandle, fullfile(FolderName, FigName));
    end
    csvwrite([FolderName,'\WL(x)_v_Width(y)__Polyfit_Surf_Intersection.csv'],[WL_roots.pos, WL_roots.neg, vec_wdt]);
    csvwrite([FolderName,'\WL(x)_v_Width(y)__Numerical_Surf_Intersection_Tolerance.csv'],[tols.x, tols.y]); 
    csvwrite([FolderName,'\WL(x)_v_Width(y)__Numerical_Surf_Intersection_MinDelta.csv'],[minskeep.WL, minskeep.widths]);
    
    
end

%% (NUMERIC-NEW) search for "root" = point at which surfaces intersect (difference is minimimal )
function [tols, minskeep, fine] = find_surface_intersections_numerical(F, SHG, titletext)

    %% USER TUNABLE PARAMETERS
    Nwidths        = 800;     % resolution of intersection curve
    Nscan_default  = 400;     % WL scan resolution
    Nscan_local    = 200;     % local scan near previous root
    WL_search_half = 25;      % continuation window [nm]
    debugPlots     = true;    % set false to disable debug figures

    %% Build interpolants (with guarded extrapolation)
    Finterp = scatteredInterpolant( ...
        F.WL(:), F.widths(:), F.vq(:), ...
        'natural', 'linear');

    SHGinterp = scatteredInterpolant( ...
        SHG.WL(:), SHG.widths(:), SHG.vq(:), ...
        'natural', 'linear');

    deltaN = @(WL,w) Finterp(WL,w) - SHGinterp(WL,w);

    %% Common domain (no blind extrapolation)
    WLrange = [ ...
        max(min(F.WL(:)),  min(SHG.WL(:))), ...
        min(max(F.WL(:)),  max(SHG.WL(:))) ];

    widthrange = [ ...
        max(min(F.widths(:)), min(SHG.widths(:))), ...
        min(max(F.widths(:)), max(SHG.widths(:))) ];

    widths_fine = linspace(widthrange(1), widthrange(2), Nwidths)';
    WL_roots    = NaN(size(widths_fine));

    %% DEBUG storage
    debug.failed_widths = [];
    debug.failed_reason = {};

    %% MAIN ROOT TRACE LOOP
    for k = 1:length(widths_fine)

        w = widths_fine(k);

        % --- Branch continuation ---
        if k > 1 && isfinite(WL_roots(k-1))
            WLtest = linspace( ...
                max(WLrange(1), WL_roots(k-1)-WL_search_half), ...
                min(WLrange(2), WL_roots(k-1)+WL_search_half), ...
                Nscan_local);
        else
            WLtest = linspace(WLrange(1), WLrange(2), Nscan_default);
        end

        % --- Evaluate deltaN safely ---
        fvals = arrayfun(@(WL) deltaN(WL,w), WLtest);

        isok  = isfinite(fvals);
        if nnz(isok) < 2
            debug.failed_widths(end+1) = w;
            debug.failed_reason{end+1} = 'insufficient finite points';
            continue
        end

        WLok = WLtest(isok);
        fok  = fvals(isok);

        s = sign(fok);
        idx = find(diff(s)~=0,1);

        if isempty(idx)
            debug.failed_widths(end+1) = w;
            debug.failed_reason{end+1} = 'no sign change';
            continue
        end

        WLbr = WLok([idx idx+1]);
        fbr  = fok([idx idx+1]);

        % --- Final guard before fzero ---
        if any(~isfinite(fbr)) || sign(fbr(1)) == sign(fbr(2))
            debug.failed_widths(end+1) = w;
            debug.failed_reason{end+1} = 'invalid bracket';
            continue
        end

        % --- High-precision root ---
        try
            WL_roots(k) = fzero(@(WL) deltaN(WL,w), WLbr);
        catch
            debug.failed_widths(end+1) = w;
            debug.failed_reason{end+1} = 'fzero failure';
        end
    end

    %% COLLECT VALID ROOTS
    ikeep = isfinite(WL_roots);

    minskeep.WL     = WL_roots(ikeep);
    minskeep.widths = widths_fine(ikeep);

    % legacy compatibility
    tols.x = minskeep.WL;
    tols.y = minskeep.widths;
    fine   = struct();

    %% MAIN RESULT PLOT
    figure;
    plot(minskeep.WL, minskeep.widths, 'k', 'LineWidth', 2);
    xlabel('WL [nm]');    ylabel('width [nm]');
    xlim([210,250]);
    title({titletext; 'High-precision numerical surface intersection'});
    grid on;

    %% DEBUG PLOTS
    if debugPlots && ~isempty(debug.failed_widths)

        % 1) Where roots failed
        figure;
        plot(minskeep.WL, minskeep.widths, 'k', 'LineWidth', 2); hold on;
        plot(WLrange(1)*ones(size(debug.failed_widths)), ...
             debug.failed_widths, 'rx');
        xlabel('WL [nm]');
        ylabel('width [nm]');
        title('Debug: widths where root finding failed');
        legend('valid root','failure location');
        grid on;

        % 2) Histogram of failure causes
        [u,~,ic] = unique(debug.failed_reason);
        counts = accumarray(ic,1);

        figure;
        bar(counts);
        set(gca,'XTickLabel',u,'XTickLabelRotation',30);
        ylabel('count');
        title('Debug: root-finding failure causes');
        grid on;
    end

end

function result = surface_intersection_advanced(F, SHG)

    %% Interpolants
    Finterp = scatteredInterpolant(F.WL(:),F.widths(:),F.vq(:),'natural','linear');
    SHGinterp = scatteredInterpolant(SHG.WL(:),SHG.widths(:),SHG.vq(:),'natural','linear');

    deltaN = @(WL,w) Finterp(WL,w) - SHGinterp(WL,w);

    widths = linspace(min(F.widths(:)),max(F.widths(:)),900)';
    WLdomain = [min(F.WL(:)), max(F.WL(:))];

    branches = {};
    branchID = 0;

    %% MULTI-BRANCH DETECTION
    for k = 1:length(widths)
        w = widths(k);
        WLtest = linspace(WLdomain(1),WLdomain(2),400);
        fvals = arrayfun(@(WL) deltaN(WL,w), WLtest);

        s = sign(fvals);
        crossings = find(diff(s)~=0);

        for c = 1:length(crossings)
            WLroot = fzero(@(WL) deltaN(WL,w), WLtest([crossings(c) crossings(c)+1]));

            if k == 1 || isempty(branches)
                branchID = branchID + 1;
                branches{branchID}.WL = WLroot;
                branches{branchID}.width = w;
            else
                % continuation: nearest branch
                d = cellfun(@(b) abs(b.WL(end)-WLroot), branches);
                [dmin, ib] = min(d);

                if dmin < 5
                    branches{ib}.WL(end+1) = WLroot;
                    branches{ib}.width(end+1) = w;
                else
                    branchID = branchID + 1;
                    branches{branchID}.WL = WLroot;
                    branches{branchID}.width = w;
                end
            end
        end
    end

    %% ERROR / CONFIDENCE ESTIMATION
    for b = 1:length(branches)
        WLb = branches{b}.WL;
        wb  = branches{b}.width;

        dfdWL = zeros(size(WLb));
        for k = 1:length(WLb)
            h = 0.05;
            dfdWL(k) = (deltaN(WLb(k)+h,wb(k)) - deltaN(WLb(k)-h,wb(k))) / (2*h);
        end

        branches{b}.sigma_WL = abs(1./dfdWL);  % sensitivity-based uncertainty
    end

    %% ANALYTIC COMPARISON (quadratic)
    for b = 1:length(branches)
        p = polyfit(branches{b}.width, branches{b}.WL, 2);
        branches{b}.analytic_fit = p;
        branches{b}.WL_fit = polyval(p, branches{b}.width);
        branches{b}.fit_error = branches{b}.WL - branches{b}.WL_fit;
    end

    %% PLOTS
    figure; hold on;
    for b = 1:length(branches)
        plot(branches{b}.WL, branches{b}.width,'LineWidth',2);
        plot(branches{b}.WL_fit, branches{b}.width,'--');
    end
    xlabel('WL [nm]');
    ylabel('width [nm]');
    title('Surface intersection branches (numeric vs analytic)');
    legend('numeric','analytic');
    grid on;

    %% OUTPUT
    result.branches = branches;
end


%% (NUMERIC OLD) search for "root" = point at which surfaces intersect (difference is minimimal )
function [tols,minskeep,fine] = find_surface_intersections_numerical_old(F,SHG,titletext) 
    deltaN_tolerance = 0.0025;

    %make refined grid 
    dWL         = 1;        %search grid refinement
    dwidth      = 1;     %search grid refinement    
    fine.vec_WLs    = [ min([F.WL(:);SHG.WL(:)]) : dWL : max([F.WL(:);SHG.WL(:)]) ]';
    fine.vec_widths = [min([F.widths(:);SHG.widths(:)]) : dwidth : max([F.widths(:);SHG.widths(:)])]';
    
    %generate *FINE* tables :  rows = width;  cols = WLs
    [fine.WLs,fine.WIDTHS] = meshgrid(fine.vec_WLs, fine.vec_widths);
    isan = find(~isnan(F.vq(:)));
    [ fine.F.vq  ] = griddata( F.WL(isan),  F.widths(isan),   F.vq(isan), fine.WLs, fine.WIDTHS );
    isan = find(~isnan(SHG.vq(:)));    
    [ fine.SHG.vq  ] = griddata( SHG.WL(isan),  SHG.widths(isan), SHG.vq(isan), fine.WLs, fine.WIDTHS );
    
    %deltaN = difference between Fundamental and SHG 
    fine.deltaN = abs(fine.F.vq - fine.SHG.vq);  %find minimum of difference
    hf2 = figure; surf( fine.vec_WLs, fine.vec_widths, fine.deltaN); xlabel('WL');ylabel('width');
    shading flat; zlabel('abs(index_1 - index_2')
    
    %[1] Method 1 -- identify curve that follows valley
    [vals,imin] = min(fine.deltaN,[],1);  isan=find(~isnan(vals)); imin = imin(isan);
    mins.x1 = fine.vec_WLs(isan); mins.y1 = fine.vec_widths(imin);
    [vals,jmin] = min(fine.deltaN,[],2);  isan=find(~isnan(vals)); jmin = jmin(isan);
    mins.x2 = fine.vec_WLs(jmin); mins.y2 = fine.vec_widths(isan);
    
    %keep only coincident points
    ikeep = find( ismember(mins.x1,mins.x2) &  ismember(mins.y1,mins.y2) );
    minskeep.WL = mins.x1(ikeep);
    minskeep.widths = mins.y1(ikeep);
    
    %[2] Method 2 -- find all XY within a tolerance DeltaN
    isel = find( fine.deltaN(:) < deltaN_tolerance );
    [rsel,csel] = ind2sub(size(fine.deltaN),isel);
    tols.x = fine.vec_WLs(csel);
    tols.y = fine.vec_widths(rsel);
    
    %Plot Phase Matching Combinations
    hf3 = figure; 
    plot(tols.x,tols.y,'*'); xlabel('WL [nm]'); ylabel('width [nm]')
    hold all; plot(mins.x1,mins.y1,'.');  
    hold all; plot(mins.x2,mins.y2,'.');  
    hold all; plot(minskeep.WL ,minskeep.widths,'-o','Linewidth',5);  xlabel('WL');
    legend(['deltaN<=',num2str(deltaN_tolerance)],'min(col)','min(row)','minSel');
    
    %polyfit -- min
    figure; plot(minskeep.WL, minskeep.widths,'.');
    [paramsA] = polyfit(minskeep.WL, minskeep.widths,2);
    Xfitrange = [min(minskeep.WL):dWL: max(minskeep.WL)]';
    hold all; plot(Xfitrange, polyval(paramsA,Xfitrange) ); 
    xlabel('WL [nm]'); ylabel('width [nm]'); 
    title({titletext;['Params: ',num2str(paramsA)]});

    %polyfit -- min
    figure; plot(tols.x, tols.y,'.');
    [paramsB] = polyfit(tols.x, tols.y, 2);
    Xfitrange2 = [min(tols.x):dWL: max(tols.x)]';
    hold all; plot(Xfitrange2, polyval(paramsB,Xfitrange2) ); 
    xlabel('WL [nm]'); ylabel('width [nm]'); 
    title({titletext;['fit to "tols" (deltaN<',num2str(deltaN_tolerance),') Params: ',num2str(paramsB)]});
    
end


%% Read File and fit surface to 3D points
function [x1q,x2q,vq,infotext] = get_surface_definition(start_path, WL, width, dN)
    %Set up / Inputs
    dirs.f  = 'G:\My Drive\Analyses(BF)\Matlab\UVQgit';
%     start_path = 'G:\My Drive\Analyses(BF)\Matlab\matDat';
 
    %Select Input Matrix
    cd(start_path)
    [filename,path] = uigetfile('*surfacedefinition*.csv','Multiselect','off');
        
    %Read and Parse Points into Lists  
    A = csvread([path,filename],1,0);
    if size(A,2)==3
        x1_list = A(:,1);  %WL
        x2_list = A(:,2);  %width
        V_list  = A(:,3);  %index    
    else
        inp.WLs     = A(1, 2:end);
        inp.widths  = A(2:end, 1);
        AA          = A(2:end,2:end);
        V_list      = reshape(AA',[],1);
        x1_list     = repmat(inp.WLs',size(AA,1),1);                %WL
        x2_list     =  sort(repmat(inp.widths',1,size(AA,2))');     %widths
        ikeep = find(x1_list~=0 & x2_list~=0 & V_list~=0);
        x1_list = x1_list(ikeep);
        x2_list = x2_list(ikeep);
        V_list = V_list(ikeep);
    end
    
    V_list = V_list+dN;    % ADD OFFSET 
    
   % figure; surf(inp.WLs,inp.widths,AA)
    
    %Interplate Surface:  Griddata
    [x1q, x2q] = meshgrid(WL, width);
    vq = griddata(x1_list', x2_list', V_list',  x1q, x2q);
    
    %plot fits v. data
    figure; 
    mesh(x1q,x2q,vq)
    hold on
    plot3(x1_list', x2_list', V_list','o');
    xlabel('WL');  ylabel('width');
    
    infotext = replace(filename(1:end-4),'_',' ');
    
    title(infotext);
    
end




%% (ANALYTIC) plot fit-surfaces 
function [WL_roots] = find_surface_intersections_analytic(F, SHG, X0, vec_WL, vec_wdt )

    %select Row (or Col) of PolyFit Params that matches X0 (reference point)
    [dummy, F.rsel]     = min(abs(F.widths(:,1) - X0.width) );  %select set of Params v. WL
    [dummy, F.csel]     = min(abs(F.WL(1,:) - X0.WL) );         %select set of Params v. WIDTH
    [dummy, SHG.rsel]   = min(abs(SHG.widths(:,1) - X0.width) );  %select set of Params v. WL
    [dummy, SHG.csel]   = min(abs(SHG.WL(1,:) - X0.WL) );         %select set of Params v. WIDTH
    
    % Fundamental -- Fit to WL  (WL Dependence)
    WL_f0.wl    = X0.WL;
    wdt_f0.wl   = X0.width;
    n_f0.wl     = F.pp_v_WL(F.rsel,3)    % 1.941420378;  %2.3133;    %fit to WL
    a1          = F.pp_v_WL(F.rsel,2)    %-0.0009252864146;     %-0.001563;     %WL dependence
    a2          = F.pp_v_WL(F.rsel,1)    %1.06E-06;            %1.07E-06;      %WL dependence ^2

    % Fundamental -- Fit to Width (WIDTH Dependence) 
    WL_f0.w     = X0.WL;
    wdt_f0.w    = X0.width;
    n_f0.w      = F.pp_v_wdt(3,F.csel)   %1.93893;     %1.69129401444789;     %fit to width
    b1          = F.pp_v_wdt(2,F.csel)   %9.57e-6;      %6.98e-4    %width dependence
    b2          = F.pp_v_wdt(1,F.csel)   %-4.92E-07;    %-4.92E-07; %width dependence ^2
    
    % SHG (TM_XX) -- Fit to WL  (WL Dependence)
    WL_s0.wl    = X0.WL;
    wdt_s0.wl   = X0.width;
    n_s0.wl     = SHG.pp_v_WL(F.rsel,3)    %1.94107;      %3.979233;    
    c1          = SHG.pp_v_WL(F.rsel,2)    %-0.003427;   %-0.01054;    
    c2          = SHG.pp_v_WL(F.rsel,1)    %= 1.33e-5;      %1.21E-05;    
    
    % SHG (TM_XX) -- Fit to Width (WIDTH Dependence) 
    WL_s0.w     = X0.WL;
    wdt_s0.w    = X0.width;
    n_s0.w      = SHG.pp_v_wdt(3,F.csel) %1.94;         %0.751;     %fit to width
    d1          = SHG.pp_v_wdt(2,F.csel) %7.26e-4;      %0.00268; 
    d2          = SHG.pp_v_wdt(1,F.csel) %-1.4e-6;      %-1.40E-06;    
    

    %plot surfaces 
    vec_n_f.wl = n_f0.wl + a1*(vec_WL-WL_f0.wl) + a2*(vec_WL-WL_f0.wl).^2;
    vec_n_s.wl = n_s0.wl + c1*(vec_WL-WL_s0.wl) + c2*(vec_WL-WL_s0.wl).^2;
    hf(2) = figure;; plot(vec_WL, [vec_n_f.wl,vec_n_s.wl]); xlabel('wavelength'); 
    title(['polyFIT: n v WL @ width=',num2str(X0.width)])
        
    vec_n_f.w = n_f0.w + b1*(vec_wdt-wdt_f0.w) + b2*(vec_wdt-wdt_f0.w).^2;
    vec_n_s.w = n_s0.w + d1*(vec_wdt-wdt_s0.w) + d2*(vec_wdt-wdt_s0.w).^2;
    figure; plot(vec_wdt, [vec_n_f.w,vec_n_s.w]); xlabel('width');
    title(['polyFIT: n v width @ WL=',num2str(X0.WL)])
        
    
    %generate Surfaces  from PolyFit parameters  
    [XX,YY] = meshgrid(vec_WL,vec_wdt);
    
    NN_f.w = n_f0.w + a2*(XX-WL_f0.w).^2 + a1*(XX-WL_f0.w) + b2*(YY-wdt_f0.w).^2 + b1*(YY-wdt_f0.w);
    NN_f.wl = n_f0.wl + a2*(XX-WL_f0.wl).^2 + a1*(XX-WL_f0.wl) + b2*(YY-wdt_f0.wl).^2 + b1*(YY-wdt_f0.wl);
    
    NN_s.w = n_s0.w + c2*(XX-WL_s0.w).^2 + c1*(XX-WL_s0.w) + d2*(YY-wdt_s0.w).^2 + d1*(YY-wdt_f0.w);
    NN_s.wl = n_s0.wl + c2*(XX-WL_s0.wl).^2 + c1*(XX-WL_s0.wl) + d2*(YY-wdt_f0.wl).^2 + d1*(YY-wdt_f0.wl);

    %plot - fundamental surface
    figure; plot3([vec_WL],[wdt_f0.wl*ones(length(vec_WL),1)],[vec_n_f.wl],'Linewidth',4); grid on;
    hold all; plot3([WL_f0.w *ones(length(vec_wdt),1)],[vec_wdt],[vec_n_f.w],'Linewidth',4); grid on;
    hold on; mesh(XX,YY,NN_f.w); hold all; mesh(XX,YY,NN_f.wl); 
    xlabel('WL'); ylabel('width'); title({'fund';['polyfit surface']});

    %plot - shg
    figure; plot3([vec_WL],[wdt_s0.wl*ones(length(vec_WL),1)],[vec_n_s.wl],'Linewidth',4); grid on;
    hold all; plot3([WL_s0.w *ones(length(vec_wdt),1)],[vec_wdt],[vec_n_s.w],'Linewidth',4); grid on;
    hold on; mesh(XX,YY,NN_s.w); hold all; mesh(XX,YY,NN_s.wl); title('SHG');
    xlabel('WL'); ylabel('width'); title({'SHG';['polyfit surface']});
    
    figure; surf(XX,YY,NN_f.w); hold all; surf(XX,YY,NN_s.w);
    xlabel('WL'); ylabel('width');

    figure; surf(XX,YY,NN_f.wl); hold all; surf(XX,YY,NN_s.wl);
    xlabel('WL'); ylabel('width');   
        
    figure; surf(XX,YY,abs(NN_f.w-NN_s.w)); 
        xlabel('WL'); ylabel('width');  title('n_f(WL,width) - n_shg(WL,width)');
    
    %find ROOT (analytically)   
    for k = 1:length(vec_wdt)
        wdt = vec_wdt(k);
        A = (a2-c2);
        B = -2*(a2*WL_f0.w - c2*WL_s0.w) + (a1 - c1);
        C = a2*WL_f0.w^2 - c2*WL_s0.w^2 + c1*WL_s0.w - a1*WL_f0.w;
        K = ( d2*(wdt-wdt_s0.w)^2 + d1*(wdt-wdt_s0.w)  ) - ...
                ( b2*(wdt-wdt_f0.w)^2 + b1*(wdt-wdt_f0.w)  );
        
        WL_roots.pos(k,1) = -B/(2*A) + (1/(2*A))*sqrt(B^2 - 4*A*(C-K) );
        WL_roots.neg(k,1) = -B/(2*A) - (1/(2*A))*sqrt(B^2 - 4*A*(C-K) );
    end
    figure; plot(WL_roots.pos,vec_wdt);  
    hold on; plot(WL_roots.neg,vec_wdt); xlim([200,500]); ylim([300,700]);
    xlabel('WL [nm]'); ylabel('width [nm]'); legend('pos.root','neg.root');

end




%% POLYFIT on each ROW
function [pp] = polyfit_manycurves_byrow(X , YY, X0, infotext)
    m_order  = 2;
    
    %X = A(1,:);         %row
    %YY = A(2:end,:);    %data - in rows
    Nrows = size(YY,1);
    
    %loop over all curves and perform Polyfit
    YYfit(:,:) = NaN(size(YY));
    for k=1:Nrows
        Y = YY(k,:);
        ifit = find(Y~=0 & ~isnan(Y));
        [pp(k,:),S] = polyfit(X(ifit)'-X0,YY(k,ifit)',m_order);
        YYfit(k,ifit) = polyval(pp(k,:),X(ifit)'-X0)';
    end
    
    %plot fits v. data
    figure; subplot(2,1,2); 
    plot(X,YY,'.'); hold on; plot(X,YYfit); xlabel(infotext)
    subplot(2,1,1); plot(X,YYfit-YY); title({'Error = Fit-Actual';['order = ',num2str(m_order)]});
   
   % outfile = [path,filename(1:end-4),'_FITparams.csv'];
   % csvwrite(outfile,pp);
end



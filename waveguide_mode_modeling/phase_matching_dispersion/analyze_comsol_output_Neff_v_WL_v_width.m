function analyze_comsol_output_Neff_v_WL_v_width()
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\helpers\waveguide');
    clear;
    close all

    dirs.d='G:\Shared drives\Corp Main\Engineering\Design & SHG Analysis (Brent)\COMSOL Modeling\Results\modal dispersion tables (for matlab)\h=335nm\';
    xfilename = '2025_12___Modal_Index_for_h=335nm.xlsx';

    grp_n_col = 10;

    S = load_comsol_sweep_cube([dirs.d,xfilename], 'Sheet10', 'a9:co70', ...
        'Format', 'column_groups', 'GroupSize', grp_n_col);

    for kg=1:S.Ngrps
        w(kg,1) = S.P2(kg);
        X =  S.P1{kg};
        Ys = S.Ys{kg};
        figure; plot(X',Ys);
        clear dYs;  dYs = zeros(size(Ys,1)-2,size(Ys,2));
        for kc=1:size(Ys,2)
            dYs(:,kc) = diff(Ys(:,kc),2);
        end
        figure; plot(X(1:end-2),dYs);
    end
    disp(S.A);



end


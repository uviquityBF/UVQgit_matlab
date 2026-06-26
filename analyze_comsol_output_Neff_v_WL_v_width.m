function analyze_comsol_output_Neff_v_WL_v_width()
    clear;
    close all

    dirs.d='G:\Shared drives\Corp Main\Engineering\Design & SHG Analysis (Brent)\COMSOL Modeling\Results\modal dispersion tables (for matlab)\h=335nm\';
    xfilename = '2025_12___Modal_Index_for_h=335nm.xlsx';

    grp_n_col = 10;

    A = xlsread([dirs.d,xfilename],'Sheet10','a9:co70');
    Ngrps = floor(size(A,2)/grp_n_col);
    cstart = mod(size(A,2),grp_n_col);
    for kg=1:Ngrps
        M{kg}=A(:,cstart+[grp_n_col*(kg-1)+1:grp_n_col*kg]);
        w(kg,1) = M{kg}(1,1);
        X =  M{kg}(2,:)';
        Ys = M{kg}(3:end,:)';
        figure; plot(X',Ys);
        clear dYs;  dYs = zeros(size(Ys,1)-2,size(Ys,2));
        for kc=1:size(Ys,2)
            dYs(:,kc) = diff(Ys(:,kc),2);
        end
        figure; plot(X(1:end-2),dYs);
    end
    disp(A);
    
    
    
end


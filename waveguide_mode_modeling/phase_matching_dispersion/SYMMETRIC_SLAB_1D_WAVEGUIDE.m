%% SYMMETRIC SLAB 1D waveguide
%
%
%

function blah = SYMMETRIC_SLAB_1D_WAVEGUIDE()
    close all;
    clear;
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\third_party');
    d           = 1e-6;    %[m] thickness of waveguide
    n1          = 2.5;      %core of waveguide
    n0          = 2.0;      %cladding (top and bottom, b/c symmetric)
    lam0_nm     = 400;      %[nm] wavelength in vacuum
    %test calculation
    modes       = get_modes_symmetric_slab(d,n1,n0,lam0_nm,1);

    %get real refractive indices
    tmp = csvread('AlN_RefractiveIndex.csv');
    aln.lam = tmp(:,1), aln.n = tmp(:,2);

    tmp = csvread('Al2O3_RefractiveIndex.csv',1,0);
    al2o3.lam = tmp(:,1)*1000, al2o3.n = tmp(:,2);
    n0 = interp1(al2o3.lam,al2o3.n,lam0_nm,'PCHIP');


    %loop over "d"  @ n_clad = 2.0
    d_list_um       = [0.4:0.01:1.6]';
    lam0_list       = [445:1:460]';
    n0              = interp1(al2o3.lam,al2o3.n,lam0_nm,'PCHIP'); % for SAPPHIRE

    for id=1:length(d_list_um)
        d = d_list_um(id,1)*1e-6;   %[um]
        legnames{id} = [num2str(d),'um'];
        for il=1:length(lam0_list)
            lam0_nm     = lam0_list(il);
            lamSHG_nm   = lam0_nm/2;
            n1          = interp1(aln.lam,aln.n,lam0_nm,'PCHIP');
            n1_SHG      = interp1(aln.lam,aln.n,lamSHG_nm,'PCHIP');
            n0          = interp1(al2o3.lam,al2o3.n,lam0_nm,'PCHIP');
            n0_SHG      = interp1(al2o3.lam,al2o3.n,lamSHG_nm,'PCHIP');

            modesFUND_sv{id,il} = get_modes_symmetric_slab(d,n1,    n0,    lam0_nm,  0);
            modesSHG__sv{id,il} = get_modes_symmetric_slab(d,n1_SHG,n0_SHG,lamSHG_nm,0);

            beta_diff_val(id,il) = mindiff_two_lists(2*modesFUND_sv{id,il}.beta_m, modesSHG__sv{id,il}.beta_m);

        end
        disp(['d progress:', num2str(id/length(d_list_um))]);
    end

    %phase mismatch v. wavelength @ many "d"
    figure; plot(lam0_list',beta_diff_val*(1e-2)'); ylabel('betaSHG[cm^-^1] - 2*betaFUND[cm^-^1]');
    xlabel('fundamental wavelength[nm]'); title('phase mismatch v. wavelength @ many slab thickness');
    legend(legnames);

    %dephasing length .... length over which you can obtain phase matching
    figure; plot(lam0_list',abs(pi./beta_diff_val')*1e3); ylabel('[mm]');
    xlabel('fundamental wavelength[nm]'); title('1pi dephasing length');
    legend(legnames);



    V1 = modesFUND_sv{id,il}.beta_m;
%     V2 = modesSHG_sv{id,il}.beta_m;
%
%     diff_val = mindiff_two_lists(V1,V2);

end

%find minimum difference between any two elements of a twos lists
function diff_val = mindiff_two_lists(V1,V2)
    if length(V1)>0 & length(V2)>0
        for k1=1:length(V1)
            [mindiffs(k1),i2] = min(abs( V2 - V1(k1) ) );
        end
        [mindiff_val,i1] = min(mindiffs);
        diff_val = V2(i2) - V1(i1);
    else
        diff_val=NaN;
    end
end


function modes = get_modes_symmetric_slab(d,n1,n0,lam0_nm,doplot)

    %conversion & calcs
    lam0            = lam0_nm*1e-9; %convert to meters
    theta_c_norm    = rad2deg(asin(n0/n1));     %critical angle rel. to boundary normal ... confine > theta_c_norm
    theta_c         = 90-rad2deg(asin(n0/n1));  %critical angle rel. to Z               ... confine < theta_x
    k0              = 2*pi/lam0;

    %number of modes
    M   = ceil( 2*d/lam0*sind(theta_c) );

    %find angles of supported modes for dielectric slab
    %(saleh,teich eq 7.2-4, p250)

    thetas_m = linspace(0,pi/2,5000)';
    RHS = sqrt(sind(theta_c)^2 ./ sin(thetas_m).^2 - 1 );
    if doplot==1    figure; plot(sin(thetas_m),RHS); ylim([0,10]);  %plot: RHS
                    ylabel('RHS, LHS'); xlabel('sin(theta_m)');
    end

    %LHS: find crossings
    xc=[];yc=[];
    for m=1:M
        LHS(:,m) = tan(pi*d/lam0*sin(thetas_m) - m*pi/2);
        i_omit = find( diff(LHS(:,m))<0 );
        LHS(i_omit,m)=NaN;
        [x0,y0,iout,jout] = intersections(sin(thetas_m),RHS, ...
                                sin(thetas_m),LHS(:,m),1);
        xc=[xc;round(x0,3)]; yc=[yc;round(y0,3)];
        if doplot==1 hold all; plot(sin(thetas_m),LHS(:,m)); %PLOT LHS
        end
    end
    xyc = unique([xc,yc],'rows');
    xc = xyc(:,1); yc=xyc(:,2);
    if doplot==1  hold on; plot(xc,yc,'o');         %PLOT solution points
        title({['Confined Modes for...'];['wavelength: ',num2str(lam0_nm),'nm;   d=',num2str(d*1e6),'um;   n_s_l_a_b=',num2str(n1),';  n_c_l_a_d=',num2str(n0)]});

    end

    %mode results
    modes.N         = length(xc);
    modes.theta_m   = asind(xc);                    %angles of supported modes
    modes.beta_m    = n1*k0*cosd(modes.theta_m);    %[1/m] propagation constants
end

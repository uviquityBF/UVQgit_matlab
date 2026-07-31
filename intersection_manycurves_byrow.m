% Intersection ofMany Curves - by Row
% Input format assumed:
%     - first row  = Yrefeerence data
%     - other rows = Y data
% ====> find intersection between each Ydata curve and Yreference curve
%
%

function intersection_manycurves_byrow()
    %Set up / Inputs
    clear; close all;
    addpath('G:\My Drive\Analyses(BF)\Matlab\UVQgit\third_party');
    dirs.f  = 'G:\My Drive\Analyses(BF)\Matlab\UVQgit';
    start_path = 'G:\My Drive\Analyses(BF)\Matlab\matDat';


    %Select Input Matrix
    cd(start_path)
    [filename,path] = uigetfile('*intersection*.csv','Multiselect','off');

    %Parse the Input Matrix
    A = csvread([path,filename]);
    Yref  = A(1,:);
    YY    = A(2:end,:);
    Ny = size(YY,1);

    %loop over all curves and perform Polyfit
    x0 = zeros(size(YY,1),1); y0=x0; iout=x0; jout=x0;
    for k=1:Ny
        Y2 = YY(k,:);
        X = [1:size(YY,2)];
        robust=1;
        [tmp1,tmp2,tmp3,tmp4] = intersections(X',Yref',X',Y2',robust);
        if size(tmp1)>0
            x0(k,1)   = tmp1;
            y0(k,1)   = tmp2;
            iout(k,1) = tmp3;
            jout(k,1) = tmp4;
        end
    end

    %plot fits v. data
    figure; plot(x0,y0,'.'); title({'Intersections'});
    %
    outdata = [x0,y0,iout,jout];
    outfile = [path,filename(1:end-4),'_INTERSECTIONS.csv'];
    csvwrite(outfile,outdata);
end

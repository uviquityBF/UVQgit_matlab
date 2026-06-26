% Polyfit to Many Curves - by Row
% Input format assumed:    
%     - first row  = X data
%     - other rows = Y data
%
%

function polyfit_manycurves_byrow()
    %Set up / Inputs
    clear; close all;
    dirs.f  = 'G:\My Drive\Analyses(BF)\Matlab\UVQgit';
    start_path = 'G:\My Drive\Analyses(BF)\Matlab\matDat';
    m_order  = 2;
    
    %Select Input Matrix
    cd(start_path)
    [filename,path] = uigetfile('*polyfit*.csv','Multiselect','off');
    
    %Parse the Input Matrix
    A = csvread([path,filename]);
    X = A(1,:);
    YY = A(2:end,:);
    Ny = size(YY,1);
    
    %loop over all curves and perform Polyfit
    YYfit(:,:) = zeros(size(YY));
    for k=1:Ny
        delta_X(k,1) = 0; %rand()*100-50;
        Y = YY(k,:);
        ifit = find(Y~=0);
        [pp(k,:),S] = polyfit(X(ifit)'-delta_X(k,1),YY(k,ifit)',m_order);
        YYfit(k,ifit) = polyval(pp(k,:),X(ifit)'-delta_X(k,1))';
    end
    
    %plot fits v. data
    figure; plot(X,YY,'.'); 
    hold on; plot(X,YYfit);
    
    figure; plot(X,YYfit-YY); title({'Error = Fit-Actual';['order = ',num2str(m_order)]});
    %
    outfile = [path,filename(1:end-4),'_FITparams.csv'];
    csvwrite(outfile,pp);
end


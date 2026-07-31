%
function test_log_linear_fitting()
    clear;
    close all;

    %% set up parameters

    a_true              = -40;                  %[dB/cm]decay rate
    Amax                =  10000;              %peak signal
    b_list              = logspace(0,4,20);     %BG counts
    NoiseMult_list      = [0.5,1,2,3,4,5];       %noise multiple
    X                   = [0:1:2000]';
    %% single case
    b=500;  NoiseMult = 4; 
	[y_true,y_fake] = generate_fake_data(X,Amax,a_true,b,NoiseMult,1);
    y_fake_log =  log10(y_fake);
    [p,S] = polyfit(X,y_fake_log,1);
    a_fit_log = p(1)*1e5
    y_fit = 10.^( polyval(p,X) );
    figure; plot(X,y_fake_log); hold all; plot(X,polyval(p,X))
    
    %% loops
    for kb=1:length(b_list)
        b = b_list(kb);
        for kn=1:length(NoiseMult_list)
            NoiseMult = NoiseMult_list(kn);

            %% get fake data and fit with Log-Linear
            [y_true,y_fake] = generate_fake_data(X,Amax,a_true,b,NoiseMult,0);
            y_fake_log =  log10(y_fake);
            [p,S] = polyfit(X,y_fake_log,1);
            a_fitlog(kb,kn)         = p(1)*1e5;
            Amax_fit_log(kb,kn)     = 10^(p(2));
            
        end
    end 
 figure; surf(NoiseMult_list , b_list/Amax   ,   a_fitlog);  ylabel('BG as fraction of True Max Signal');
 xlabel('Noise Multiple (1= 1 sigma)');   zlabel('Log-Linear Fit Estimate');
 title({'Fit Estimateion by Log-linear fit,';['based on generated data with true value of ',num2str(a_true)]});
 figure; surf(Amax_fit_log);
    
end  %MAIN

%% generate fake data
function [y_true,y_fake] = generate_fake_data(X,Amax,a_true,b,NoiseMult,do_plot);
    Ns                  = size(X,1);
    y_true              = Amax*10.^((a_true/10)*X*1e-4);
    y_fake_noisefree    = (y_true + b)  
    y_fake              = y_fake_noisefree + randn(Ns,1).*sqrt(y_fake_noisefree)*NoiseMult;
    if do_plot==1
        figure; plot(X,[y_fake,y_fake_noisefree]);  %ylim([0,Amax]);
    end
end





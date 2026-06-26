%% This code is to test SHG efficiency in AlN linear WG
%% Study the scattered light collection
%% phase matching bandwidth
%% Checked on 11/17/2025

clear all

c_const = 2.9979e8; % vacuum light speed (m/s)
eps_0 = 8.854e-12; % vacuum permittivity in m^(-3)kg^(-1)s^4*A^2

P_SH_out =[];
P_SH_out2 =[];

wl1 = 0.446e-6;  % optical pump FM wavelength in m
wl1_0 = 0.446e-6;
for index = 1:length(wl1)
freq1 = c_const/wl1(index); % optical FM freq in Hz
omega1 = 2*pi*freq1; % optical FM angular freq in rad


wl2 = wl1(index)/2; %0.23e-6;  % generated SH wavelength in m
freq2 = c_const/wl2; %  SH freq in Hz
omega2 = 2*pi*freq2; % optical SH angular freq in rad


n1 = 2.04;
ng1 = 2.3; % group index for FM mode
k1 = omega1*n1/c_const;

n2 = n1; % effective index for SH mode
ng2 = 4.76; % group index for SH mode
k2 = omega2*n2/c_const;


dbeta = 1; % phase matched (use 1 instead of 0 to avoid singularity)

loss1_db = 10e2; % WG loss at FM 5 dB/cm = 500 dB/m
alpha1 = loss1_db/4.343; % change to loss per m

loss2_db = 10*loss1_db; % WG loss at SH
alpha2 = loss2_db/4.343; % change to loss per m

dalpha = alpha2/2-alpha1;

L_SHG = 2.0e-3; % length of SHG WG in m
dL0 = 0.025e-3;
L = 0.05e-3:dL0:L_SHG; % SHG WG length in m

ol_factor = 0.048e6; % 1/m  % nonlinear mode overlap factor: 0.048e6 for TM0&TM04 (1/um)--W300nm, 0.037e6---W600 nm
deff = 7.0e-12; % nonlinear coefficient (m/V) for AlN deff = d33 = 4.7 - 11 pm/V

yeta = 2*omega1^2/(eps_0*c_const^3*n1^2*n2)*deff^2*ol_factor^2; %normalized conversion efficiency W*m^(-2)
yeta2 = yeta*1e-4; % convert to W*cm^(-2)

g = sqrt(yeta);

PM_factor = exp(-alpha2*L).*((exp(dalpha*L)-1).^2+4*exp(dalpha*L).*(sin(dbeta*L/2)).^2)./((dalpha*L).^2+(dbeta*L).^2);
PM_factor2 = (sin(dbeta*L/2)./(dbeta*L/2)).^2; % ideal lossless case

P_factor = PM_factor.*L.^2; % Power-independent L dependent factor

CE = yeta*L.^2.*PM_factor;  % SHG efficiency 1/W
P_FM = 200e-3; % pump power in W
Input_CE = 0.6; % input coupling efficiency
P_FM = P_FM*Input_CE;
P_SH = CE*P_FM.^2*1e6;  % SH power in uW

%%solve ode45 for SHG coupled mode equations (CMEs)
j = sqrt(-1);
fun = @(t, x)[j*g*conj(x(1))*x(2)*exp(j*dbeta*t)-alpha1/2*x(1);       % 1st ODE
              j*g*x(1)^2*exp(-j*dbeta*t)-alpha2/2*x(2)];    % 2nd ODE


zspan = [0 L_SHG];

y0 = sqrt(P_FM);                           % initial value of y(t)
z0 = 0;                      % initial value of z(t)
x0 = [y0; z0];
[Z, B] = ode45(fun, zspan, x0);     % Solver assigns solutions to Z array
b1 = B(:,1);                         % numerical solution for y(t)
b2 = B(:,2);                         % numerical solution for z(t)


P_pump_out = abs(b1).^2*1e6;
P_SH_out = P_SH; % analytical solution
P_SH_out2 = abs(b2).^2*1e6; % CME solution

%% Solve for Scattered SHG Power 
% Ps =zeros(1,length(L));
% Pss = Ps;
% 
% for i = 2:length(L)
% 
%     Ps(i) = Ps(i-1)+ P_SH(i-1)*(1-exp(-alpha2*dL));
%     Pss(i) = P_SH(i-1)*(1-exp(-alpha2*dL));
% end

% P_s_out(index) = Ps(end);

% P_total_out(index) = P_SH_out(index) + P_s_out(index);


end

figure(1)
plot(L*1e3,P_SH_out,'r') % plot analytic SHG power
hold on
plot(Z*1e3,P_SH_out2,'b') % plot CME SHG power
xlabel('WG length (mm)')
ylabel('SHG power (uW)')




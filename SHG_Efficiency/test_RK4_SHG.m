% Gemini Code 
% for testing the RK4 Engine

function test_RK4_SHG()
    clear; close all;

    % --- 1. CONFIGURATION (Standard Case) ---
    c_c0 = 3e8; 
    dz = 0.25;               % [um] Step size
    L_total = 2000;          % [um] Waveguide Length
    z = 0:dz:L_total;
    
    % Material/Mode Parameters (TM04 example)
    Pp = 0.1;                % [W] Input Pump Power
    g = 1e-4;                % [um^-1] Nonlinear coupling (sqrt of prefactor)
    dk_constant = 0.005;     % [um^-1] Phase mismatch (Delta k)
    
    % Loss Parameters (dB/cm to um^-1 conversion)
    loss_p_dBcm = 500;         % Pump loss
    loss_s_scat_dBcm = 100;  % SHG scattering (our signal)
    loss_s_abs_dBcm = 20;    % SHG absorption (dead end)

    a0 = (loss_p_dBcm / 4.3429) * 1e-4;
    a3_scat = (loss_s_scat_dBcm / 4.3429) * 1e-4;
    a3_abs = (loss_s_abs_dBcm / 4.3429) * 1e-4;
    a3_total = a3_scat + a3_abs;

    
    
    
    
    
    % --- 2. INITIALIZATION ---
    % State Vector Y = [A1_real; A1_imag; A3_real; A3_imag; P_scat; Int_dk]
    % We separate real/imag to keep the RK4 math straightforward
    Y = zeros(6, length(z));
    Y(:,1) = [sqrt(Pp); 0; 0; 0; 0; 0]; 

    % --- 3. RK4 INTEGRATION ---
    for i = 1:(length(z)-1)
        k1 = shg_derivatives(z(i),        Y(:,i),      a0, a3_total, a3_scat, g, dk_constant);
        k2 = shg_derivatives(z(i)+dz/2,   Y(:,i)+k1*dz/2, a0, a3_total, a3_scat, g, dk_constant);
        k3 = shg_derivatives(z(i)+dz/2,   Y(:,i)+k2*dz/2, a0, a3_total, a3_scat, g, dk_constant);
        k4 = shg_derivatives(z(i)+dz,     Y(:,i)+k3*dz,   a0, a3_total, a3_scat, g, dk_constant);
        
        Y(:,i+1) = Y(:,i) + (dz/6)*(k1 + 2*k2 + 2*k3 + k4);
    end

    % --- 4. POST-PROCESSING ---
    P1 = Y(1,:).^2 + Y(2,:).^2;
    P3_guided = Y(3,:).^2 + Y(4,:).^2;
    P3_scat = Y(5,:); % Already in Watts
    
    % --- 5. PLOTTING ---
    figure('Color','w','Position',[100 100 800 600]);
    subplot(2,1,1);
    semilogy(z, P1*1e6, 'k', 'LineWidth', 1.5); hold on;
    semilogy(z, P3_guided*1e6, 'b', 'LineWidth', 1.5);
    semilogy(z, P3_scat*1e6, 'r--', 'LineWidth', 1.5);
    grid on; ylabel('Power [\muW]'); legend('Pump','SHG Guided','SHG Scattered');
    title('RK4 Solver: Absolute Power Channels');

    subplot(2,1,2);
    plot(z, P3_scat ./ (P3_guided + 1e-12), 'm', 'LineWidth', 1.5);
    grid on; ylabel('Ratio'); xlabel('z [\mum]');
    title('Ratio: Scattered / Guided SHG');
end

function dY = shg_derivatives(~, Y, a0, a3_total, a3_scat, g, dk)
    % Unpack State
    A1 = Y(1) + 1i*Y(2);
    A3 = Y(3) + 1i*Y(4);
    accum_phase = Y(6); % Integrated Delta k
    
    % 1. Pump Derivative
    dA1 = -0.5 * a0 * A1;
    
    % 2. SHG Derivative (The source term uses the accumulated phase)
    % dA3/dz = -alpha/2 * A3 + i * g * A1^2 * exp(-i * phase)
    dA3 = -0.5 * a3_total * A3 + 1i * g * (A1^2) * exp(-1i * accum_phase);
    
    % 3. Scattered Power Derivative (Power Bucket)
    % dPscat/dz = alpha_scat * |A3|^2
    dPscat = a3_scat * (abs(A3)^2);
    
    % 4. Phase Accumulation Derivative
    dPhase = dk; % In future, dk can be dk(z) here
    
    % Pack Derivative
    dY = [real(dA1); imag(dA1); real(dA3); imag(dA3); dPscat; dPhase];
end
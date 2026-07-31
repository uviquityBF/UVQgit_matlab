function [DeRate_guided, DeRate_scatt] = shg_analytic_derate(dk_per_cm, L_um, a0_dBcm, a3_dBcm)
% Analytic (closed-form) derate factors for phase-mismatched, lossy SHG --
% the fraction of the ideal (zero-detuning, zero-loss) conversion efficiency
% that survives, for guided and scattered SH power respectively.
%
% Salvaged from deprecated_scripts\shg_efficiency_2024-25\SHG_Efficiency_Model_linear_analyticOnly.m
% (function PhaseMismatch_and_Loss_Derate, "Wei's formula" / DeRate_Wei + DeRate_Scatt).
% NOT currently wired into shg_efficiency\SHG_Design_Suite.m, which uses only
% numeric RK4 integration -- that analytic path was computed-but-unused as of
% SHG_Design_Suite_v2_2 and dropped entirely from v3.0 onward. Kept here as a
% ready-to-use reference/cross-check per an earlier request to preserve it,
% not auto-restored into the active script (that's a feature decision, not a
% cleanup one).
%
%   dk_per_cm : phase mismatch [1/cm] = 4*pi/lam1_cm * dn * delta_lambda
%   L_um      : waveguide length [um]
%   a0_dBcm   : pump loss [dB/cm]
%   a3_dBcm   : SH loss [dB/cm] (combined guided+scattered loss, as in the source model)

    a0_cm = a0_dBcm / 4.3429;   % dB/cm -> Np/cm
    a3_cm = a3_dBcm / 4.3429;
    L_cm  = L_um * 1e-4;

    dkL = dk_per_cm * L_cm;
    a3L = a3_cm * L_cm;
    a0L = a0_cm * L_cm;
    daL = a3L/2 - a0L;

    if dkL == 0
        if daL == 0
            DeRate_guided = exp(-a3L) * sinc(dkL/2/pi)^2;
        else
            DeRate_guided = exp(-a3L) * ((exp(daL)-1)^2 + 4*exp(daL)*(sin(dkL/2))^2) / (dkL^2 + daL^2);
        end
    else
        DeRate_guided = exp(-a3L) * ((exp(daL)-1)^2 + 4*exp(daL)*(sin(dkL/2))^2) / (dkL^2 + daL^2);
    end

    DeRate_scatt = ((0.5*a3L)^2) * ((exp(-a0L)-1)^2 + 4*exp(-a0L)*(sin(dkL/2))^2) ...
                   / ((dkL^2 + daL*a0L)^2 + ((daL*dkL) - (a0L*dkL))^2);
end

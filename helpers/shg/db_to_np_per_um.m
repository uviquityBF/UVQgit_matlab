function a_np_per_um = db_to_np_per_um(a_dBcm)
% Converts a loss coefficient from dB/cm to Np/um (field-amplitude decay rate).
% For a dB/cm/mm gradient, multiply this function's result by 1e-3 (an extra
% length-unit factor for the "per mm" of the gradient's denominator).
    a_np_per_um = (a_dBcm / 4.3429) * 1e-4;
end

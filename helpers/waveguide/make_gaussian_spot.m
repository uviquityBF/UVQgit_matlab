function spot = make_gaussian_spot(cfg)
% Coherent elliptical Gaussian input spot defined by its mode field diameter (MFD).
%
% USAGE:
%   spot = make_gaussian_spot(cfg);       % uses cfg.gauss (see below)
%
% cfg.gauss fields
%   mfd_x_um, mfd_y_um   mode field diameter along x / y (um) = 1/e^2 INTENSITY
%                        diameter (2*w, w = 1/e^2 intensity radius = 1/e amplitude
%                        radius). Equivalent FWHM(intensity) = MFD*sqrt(ln2/2) = 0.589*MFD.
%   m2                   beam quality (default 1). A single coherent field cannot have
%                        M2 > 1, so M2 only changes the divergence: it scales the Rayleigh
%                        range (z_R = pi*w0^2/(M2*lambda)) when the beam is defocused.
%   lambda_um            vacuum wavelength, used only for defocus (default 0.45)
%   window_um, n_pts     native-grid window (um, default 12) and points per axis (default 1024)
% The spot is at its waist (flat phase) on the grid; defocus / offset are applied
% later by modal_input_field (cfg.defocus_um, cfg.beam_offset_um).
%
% Returns the same struct as load_input_spot:
%   zX, zY (um), E_native, P_in_total, cx_in, cy_in, lambda_um (propagation
%   wavelength, already multiplied by M2), n_medium, name

    g = struct();  if isfield(cfg, 'gauss'), g = cfg.gauss; end
    mfd_x  = getopt(g, 'mfd_x_um', 1.0);
    mfd_y  = getopt(g, 'mfd_y_um', mfd_x);
    m2     = getopt(g, 'm2', 1);
    lam    = getopt(g, 'lambda_um', 0.45);
    W      = getopt(g, 'window_um', 12);
    n      = getopt(g, 'n_pts', 1024);
    if max(mfd_x, mfd_y) > W / 3
        warning('make_gaussian_spot:window', 'MFD (%.2f um) is large compared with the %.1f um native window.', ...
                max(mfd_x, mfd_y), W);
    end

    x = linspace(-W / 2, W / 2, n);   y = x;
    [X, Y] = meshgrid(x, y);
    E = exp(-4 * (X.^2 / mfd_x^2 + Y.^2 / mfd_y^2));       % amplitude; intensity = exp(-8 r^2/MFD^2)

    spot.zX = x;  spot.zY = y;  spot.E_native = E;
    spot.P_in_total = sum(abs(E(:)).^2) * (x(2) - x(1)) * (y(2) - y(1));
    spot.cx_in = 0;  spot.cy_in = 0;
    spot.lambda_um = lam * m2;          % effective propagation wavelength (M2 scales the divergence)
    spot.n_medium = 1;
    spot.name = sprintf('Gaussian MFD %.2f x %.2f um (M2=%.2f)', mfd_x, mfd_y, m2);
end

function v = getopt(s, name, default)
    if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

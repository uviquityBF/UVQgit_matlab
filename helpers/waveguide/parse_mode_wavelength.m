function lam_um = parse_mode_wavelength(name)
% Wavelength (um) that a COMSOL mode file was solved at, read from the "wl=<number>nm" tag in
% its name, e.g.  modeData_wl=450nm__N=20__h=350nm__w=2500nm.csv  ->  0.450.
% The path is ignored. Returns NaN if there is no tag.
    [~, nm, ext] = fileparts(name);
    t = regexp([nm ext], 'wl=(\d+\.?\d*)nm', 'tokens', 'once');
    if isempty(t), lam_um = NaN; else, lam_um = str2double(t{1}) / 1000; end
end

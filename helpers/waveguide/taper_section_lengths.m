function [L_sec, z_edges] = taper_section_lengths(widths_um, profile_fn, L_um)
% Length of each staircase section of a taper, for sections centred on the
% given (exported) widths. widths_um must be sorted from widest to narrowest.
%
% USAGE:
%   [L_sec, z_edges] = taper_section_lengths(widths_um, profile_fn, L_um)
%
% profile_fn(s) maps normalized taper position s in [0,1] to the REMAINING
% fraction of the width change: f(0) = 1 (start width), f(1) = 0 (end width),
% monotonically decreasing. Examples:
%   linear                  @(s) 1 - s
%   fast narrowing first    @(s) (1 - s).^2
%   slow narrowing first    @(s) 1 - s.^2
%   S-curve (cosine)        @(s) 0.5*(1 + cos(pi*s))
%
% Section i covers the z-range between the boundaries placed at the mid-widths
% of neighbouring exported widths. L_sec sums to L_um; z_edges has n+1 entries.

    w = widths_um(:);  n = numel(w);
    if any(diff(w) >= 0), error('widths_um must be strictly decreasing.'); end
    if n == 1, L_sec = L_um; z_edges = [0 L_um]; return; end

    sg = linspace(0, 1, 20001)';
    fg = profile_fn(sg);  fg = fg(:);
    if abs(fg(1) - 1) > 1e-6 || abs(fg(end)) > 1e-6
        error('profile_fn must satisfy f(0)=1 and f(1)=0 (got %.3g, %.3g).', fg(1), fg(end));
    end
    if any(diff(fg) > 1e-9), error('profile_fn must be monotonically decreasing.'); end

    frac = (w - w(end)) / (w(1) - w(end));              % 1 ... 0
    mid  = (frac(1:end - 1) + frac(2:end)) / 2;         % boundary widths (as fractions)
    [fu, ia] = unique(fg);                              % ascending f for interpolation
    s_mid = interp1(fu, sg(ia), mid, 'linear');         % position where the taper hits each mid-width
    s_edges = [0; s_mid(:); 1];
    z_edges = L_um * s_edges';
    L_sec = diff(z_edges);
end

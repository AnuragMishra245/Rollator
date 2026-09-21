function fig = stability_balance_map(src, opts)
%STABILITY_BALANCE_MAP  Thesis figure: rearward/forward tipping-capacity
% balance over X_r and X_cm at a fixed wheelbase L. Interpolated straight
% from the per-design Psi_r/Psi values stability_grid_pipeline.m already
% computed, not recomputed, so it can't disagree with the workbook.
%
% Colour = Psi_r/Psi. 0.5 = balanced (the goal), >0.5 forward tipping binds,
% <0.5 backward tipping binds. The 0.5 contour is the balance line
% X_cm*(X_r) = (X_r + L/2)(1 + Fz,max/(2 m_r g)); Psi is constant over the
% whole plane (only depends on L). R26, if given, marks the heel-clearance
% limit; past it the data's still real, just not adoptable.
%
%   stability_balance_map                              % reads stability_grid_out.xlsx
%   stability_balance_map(T)                           % from the pipeline table
%   stability_balance_map("stability_grid_out.xlsx", 'Save',"balance_map.png")
%   stability_balance_map(..., 'ShowTitle',false, 'Save',"balance_map.pdf")   % thesis form
%
% OPTS
%   .Save          one or more paths to save (e.g. ["fig.svg","fig.png"]);
%                  vector formats (.svg/.pdf/.eps) export as vector, others
%                  (.png/.jpg) export at 300 dpi                            [""]
%   .Fraction      "fraction" | "percent"   colour-axis units          ["percent"]
%   .FlipColor     swap the blue / orange sides                          [false]
%   .ShowTitle     draw the title (turn OFF and use a LaTeX caption)      [true]
%   .Font          axis font                                 ["Times New Roman"]
%   .L .mr .Fzmax .g .a_max   pinned parameters (overridden by columns of src)
%   .Xr_lim .Xcm_lim   plot window, m            [-1.20 -0.85] / [-0.90 -0.20]
%   .R26           heel-clearance limit to mark, if any (m)              [NaN]

arguments
    src = ""
    opts.Save          (1,:) string  = ""
    opts.Fraction      (1,1) string  = "percent"
    opts.FlipColor     (1,1) logical = false
    opts.ShowTitle     (1,1) logical = true
    opts.Font          (1,1) string  = "Times New Roman"
    opts.L             (1,1) double  = 1.20
    opts.mr            (1,1) double  = 32.964
    opts.a_max         (1,1) double  = 1.5
    opts.Fzmax         (1,1) double  = 300
    opts.g             (1,1) double  = 9.81
    opts.Xr_lim        (1,2) double  = [-1.20 -0.85]
    opts.Xcm_lim       (1,2) double  = [-0.90 -0.20]
    opts.R26           (1,1) double  = NaN   % heel-clearance limit (eq. 4.2), if any
end

% load the pipeline table
if istable(src)
    T = src;
elseif strlength(string(src)) > 0 && isfile(char(src))
    T = readtable(char(src), 'TextType', 'string');
else
    def = fullfile(fileparts(mfilename('fullpath')), 'grid_results', 'stability_grid_out.xlsx');
    assert(isfile(def), 'stability_balance_map:noData', ...
        'No table given and no workbook found at %s', def);
    T = readtable(def, 'TextType', 'string');
end
req = ["Xr","Xcm","Psi","Psi_r"];
for f = req
    assert(ismember(f, string(T.Properties.VariableNames)), ...
        'stability_balance_map:missingCol', 'input table is missing column "%s".', f);
end
for f = ["L","mr","a_max","Fzmax","g"]
    opts.(f) = pickcol(T, f, opts.(f));
end

L = opts.L;  mr = opts.mr;  amax = opts.a_max;  Fzm = opts.Fzmax;  g = opts.g;
kbal = 1 + Fzm / (2*mr*g);
R26  = opts.R26;
tf   = char(opts.Font);

pct = strcmpi(opts.Fraction, "percent");
sc  = 1;  ctr = 0.5;  cbl = 'rearward share of tipping capacity,   \Psi_r / \Psi';
if pct
    sc = 100;  ctr = 50;  cbl = 'rearward share of tipping capacity,   \Psi_r / \Psi   (%)';
end

% interpolate the real per-design data onto a fine mesh
FI = scatteredInterpolant(T.Xr, T.Xcm, T.Psi_r ./ T.Psi, 'linear', 'none');

xr  = linspace(opts.Xr_lim(1), opts.Xr_lim(2), 600);
xcm = linspace(opts.Xcm_lim(1), opts.Xcm_lim(2), 600);
[XR, XCM] = ndgrid(xr, xcm);
Xf = XR + L;

FRp = FI(XR, XCM) * sc;                             % NaN outside the data's convex hull already
feas = (XCM > XR) & (XCM < Xf);                     % support polygon: the only hard mask
FRp(~feas) = NaN;

lo = min(FRp(:), [], 'omitnan') - 0.015*sc;
hi = max(FRp(:), [], 'omitnan') + 0.015*sc;

% figure
fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [90 50 900 900]);
if isprop(fig, 'Theme'); try; fig.Theme = 'light'; catch; end; end
ax = axes(fig);  hold(ax, 'on');  ax.SortMethod = 'childorder';

contourf(ax, XR, XCM, FRp, linspace(lo, hi, 40), 'LineColor', 'none');
colormap(ax, diverge_twoslope(256, lo, ctr, hi, opts.FlipColor));
clim(ax, [lo hi]);

step = 0.10*sc;
lv = (ceil(lo/step)*step : step : floor(hi/step)*step);
lv = lv(abs(lv - ctr) > 1e-9);
[cc, hcl] = contour(ax, XR, XCM, FRp, lv, 'LineColor', [0.38 0.38 0.38], 'LineWidth', 0.45);
clabel(cc, hcl, 'FontName', tf, 'FontSize', 8, 'Color', [0.33 0.33 0.33], 'LabelSpacing', 460);

% R26 heel-clearance limit, if given and inside the window
if ~isnan(R26) && R26 > opts.Xr_lim(1) + 1e-9 && R26 < opts.Xr_lim(2) - 1e-9
    xline(ax, R26, 'k:', 'LineWidth', 1.4);
    % placed well below the balance line's crossing of R26, into the orange
    % (Psi_f > Psi_r) region, so it doesn't compete with the red label
    yLabel = opts.Xcm_lim(2) - 0.60*(opts.Xcm_lim(2) - opts.Xcm_lim(1));
    text(ax, R26 - 0.014, yLabel, sprintf('X_r = %.2f m  (heel-clearance limit)', R26), ...
        'FontName', tf, 'FontSize', 8, 'Rotation', 90, 'BackgroundColor', [1 1 1], ...
        'Margin', 1, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
else
    R26 = opts.Xr_lim(2);   % nothing to mark -> the balance line just runs the full window
end

xrl = linspace(opts.Xr_lim(1), opts.Xr_lim(2), 600);
xcl = kbal * (xrl + L/2);
inwin = xcl >= opts.Xcm_lim(1) & xcl <= opts.Xcm_lim(2);
segf  = inwin & xcl > xrl & xcl < xrl + L & xrl <= R26;
plot(ax, xrl(inwin), xcl(inwin), '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.0);
plot(ax, xrl(segf),  xcl(segf),  '-',  'Color', [0.78 0.05 0.09], 'LineWidth', 3.6);

% aspect + limits (equal, so the balance line reads at its true slope)
daspect(ax, [1 1 1]);
axis(ax, [opts.Xr_lim opts.Xcm_lim]);

% annotations
ang = atan2d(kbal, 1);
nrm = hypot(1, kbal);
xrS = xrl(segf);
if ~isempty(xrS)
    xm0 = xrS(1) + 0.46*(xrS(end) - xrS(1));
    xa  = xm0 - 0.045*kbal/nrm;  ya = kbal*(xm0 + L/2) + 0.045/nrm;   % offset into the blue side
    text(ax, xa, ya, '\Psi_r = \Psi_f', 'Rotation', ang, 'FontName', tf, ...
        'FontSize', 11, 'Color', [0.78 0.05 0.09], 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

text(ax, 0.035, 0.975, {'\Psi_r > \Psi_f', 'forward-tipping limited'}, 'Units', 'normalized', ...
    'FontName', tf, 'FontSize', 9.5, 'Color', [0.10 0.22 0.42], 'VerticalAlignment', 'top');
text(ax, 0.965, 0.03, {'\Psi_f > \Psi_r', 'backward-tipping limited'}, 'Units', 'normalized', ...
    'FontName', tf, 'FontSize', 9.5, 'Color', [0.50 0.22 0.06], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% axes cosmetics
box(ax, 'on');  grid(ax, 'on');  ax.Layer = 'top';
ax.GridColor = [0.5 0.5 0.5];  ax.GridAlpha = 0.10;
ax.FontName = tf;  ax.FontSize = 11;  ax.LineWidth = 0.8;
ax.XColor = [0 0 0];  ax.YColor = [0 0 0];
tstep = 0.05;
ax.XTick = ceil(opts.Xr_lim(1)/tstep)*tstep  : tstep : floor(opts.Xr_lim(2)/tstep)*tstep;
ax.YTick = ceil(opts.Xcm_lim(1)/tstep)*tstep : tstep : floor(opts.Xcm_lim(2)/tstep)*tstep;
xlabel(ax, 'rear-axle position   X_r   (m)', 'FontName', tf, 'FontSize', 12);
ylabel(ax, 'CoM fore-aft position   X_{cm}   (m)', 'FontName', tf, 'FontSize', 12);
if opts.ShowTitle
    title(ax, {'Rearward / forward tipping-capacity balance', ...
        sprintf('L = %.3f m fixed,  a_{max} = %.2g m/s^2', L, amax)}, ...
        'FontName', tf, 'FontSize', 12, 'FontWeight', 'bold');
end

cb = colorbar(ax);
cb.Label.String = cbl;  cb.Label.FontName = tf;  cb.Label.FontSize = 11;
cb.FontName = tf;  cb.FontSize = 9.5;  cb.LineWidth = 0.8;
t10 = 0.1*sc;
cb.Ticks = unique(sort([ (ceil(lo/t10)*t10 : t10 : floor(hi/t10)*t10), ctr ]));

drawnow;

for i = 1:numel(opts.Save)
    p = opts.Save(i);
    if strlength(p) == 0; continue; end
    [~, ~, ext] = fileparts(char(p));
    if any(strcmpi(ext, {'.svg','.pdf','.eps'}))
        exportgraphics(fig, p, 'ContentType', 'vector');
    else
        exportgraphics(fig, p, 'Resolution', 300);
    end
end

end % stability_balance_map



function v = pickcol(T, name, dflt)
    if ~isempty(T) && ismember(name, string(T.Properties.VariableNames))
        v = T.(name)(1);
    else
        v = dflt;
    end
end


function cm = diverge_twoslope(n, lo, mid, hi, flip)
%DIVERGE_TWOSLOPE  orange (low) -> off-white (mid) -> blue (high), with the
% white pinned at `mid` even when [lo, hi] is not symmetric about it.
    if nargin < 5; flip = false; end
    loC = [0.800 0.350 0.110];   % Psi_r/Psi low  : forward-capacity dominates
    midC= [0.972 0.972 0.965];
    hiC = [0.130 0.400 0.670];   % Psi_r/Psi high : rearward-capacity dominates
    if flip; tmp = loC; loC = hiC; hiC = tmp; end
    f  = (mid - lo) / (hi - lo);
    f  = min(max(f, 0.06), 0.94);
    nlo = max(round(n*f), 2);  nhi = n - nlo;
    c1 = [linspace(loC(1),midC(1),nlo).' linspace(loC(2),midC(2),nlo).' linspace(loC(3),midC(3),nlo).'];
    c2 = [linspace(midC(1),hiC(1),nhi).' linspace(midC(2),hiC(2),nhi).' linspace(midC(3),hiC(3),nhi).'];
    cm = [c1; c2];
end

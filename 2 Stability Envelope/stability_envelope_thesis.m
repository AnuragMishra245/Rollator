function result = stability_envelope_thesis(design, opts)
%STABILITY_ENVELOPE_THESIS  Dynamic tipping envelope, trapezoid check, area Psi.
% Implements Bachelorarbeit.pdf Sec. 3.3.2/3.3.3 (tipping boundary, envelope)
% and Sec. 3.4.1/3.4.2 (Fz,c, non-degeneracy, area Psi). The equations are
% below, each one tagged with its thesis equation number.
%
% Press Run (or call with no args) to plot one design via the manual input
% section near the bottom. Called from the pipeline with Plot=false, only
% the core computation runs. Programmatic use:
%   d = struct('mr',19.624,'H',1.075,'Xr',-0.300,'Xcm',0.001672,'Xf',0.300,'Zcm',0.211226);
%   r = stability_envelope_thesis(d, 'Name',"My concept", 'a_max',1.5);
%
% Sign convention (Ackermann & Mombaur Fig. 1(d) / thesis Fig. 3.9): origin O
% is the ground point below the handles; Xr<0<Xf. F_x>0 = toward the user
% (backward-tipping risk). F_z>0 = downward. a_x>0 = forward. theta>0 = front
% wheel downhill (theta_deg in degrees, everything else SI).
%
% Not implemented here: shape coordinates (Sec. 3.3.5, "3 Parameter Study"),
% the R/B split + ellipse-to-boundary margin (Sec. 3.4.3-3.4.5, "4 Stability
% Metric"), the 3-stage optimization (Sec. 3.4.6, "5 Optimization"). Those
% files reuse result.fn (the four boundary functions) instead of re-deriving.

arguments
    design = []                            % [] runs the manual input section below

    opts.theta_deg   (1,1) double = 0       % ground incline, degrees (front downhill +)
    opts.a_max       (1,1) double = 1.5     % prescribed |a_x| bound, m/s^2
    opts.Fzmax       (1,1) double = 300     % fill top edge AND Psi upper limit, N
    opts.FzFillMin   (1,1) double = 0       % fill bottom edge AND Psi lower limit, N (thesis fixes this at 0)
    opts.g           (1,1) double = 9.81
    opts.crr         (1,1) double = 0.05    % pass-through only, crr doesn't enter Eqs 3.28/3.29/3.48/3.49/3.51,
                                            % kept just so the parameter table has all columns
    opts.Plot        (1,1) logical = true
    opts.Visible     (1,1) logical = true
    opts.RunAsserts  (1,1) logical = true   % accepted but unused, keeps old pipeline calls with 'RunAsserts',false working
    opts.Verbose     (1,1) logical = true
    opts.Name        (1,1) string  = "design"
end

if isempty(design)
    manualInputSection();
    result = [];
    return
end

result = coreRun(design, opts, opts.Name);

end % stability_envelope_thesis


%% Equations, one-to-one with Bachelorarbeit.pdf Ch. 3
% Pure functions only, no plotting/validation/pipeline glue. Each formula is
% copied directly from the cited equation.

% Eq. (3.28), backward (rear) tipping boundary: F_z as a function of F_x,
% general incline theta, prescribed acceleration a_x.
function Fz = fzBack(Fx, ax, d, g, theta)
    Fz = ( d.mr*g*cos(theta)*(d.Xcm - d.Xr) + d.mr*g*d.Zcm*sin(theta) ...
           - d.H.*Fx - d.mr*d.Zcm.*ax ) / d.Xr;
end

% Eq. (3.29), forward (front) tipping boundary.
function Fz = fzForw(Fx, ax, d, g, theta)
    Fz = ( d.mr*g*cos(theta)*(d.Xcm - d.Xf) + d.mr*g*d.Zcm*sin(theta) ...
           - d.H.*Fx - d.mr*d.Zcm.*ax ) / d.Xf;
end

% Eq. (3.48), Eq. (3.28) solved for F_x (the plotting form; also used
% directly as F+, the rear-tipping edge at a_x=+a_max, thesis Eq. 3.41).
function Fx = fxBack(Fz, ax, d, g, theta)
    Fx = ( d.mr*g*cos(theta)*(d.Xcm - d.Xr) + d.mr*g*d.Zcm*sin(theta) ...
           - d.mr*d.Zcm.*ax - d.Xr.*Fz ) / d.H;
end

% Eq. (3.49), Eq. (3.29) solved for F_x (also F-, thesis Eq. 3.42, at a_x=-a_max).
function Fx = fxForw(Fz, ax, d, g, theta)
    Fx = ( d.mr*g*cos(theta)*(d.Xcm - d.Xf) + d.mr*g*d.Zcm*sin(theta) ...
           - d.mr*d.Zcm.*ax - d.Xf.*Fz ) / d.H;
end

% Eq. (3.43), crossing point F_z,c of the two a_x=+-a_max boundaries.
function Fzc = crossingPointFzc(d, g, theta, amax)
    Fzc = -d.mr*g*cos(theta) + 2*d.mr*d.Zcm*amax / (d.Xf - d.Xr);
end

% Eq. (3.45), non-degeneracy (trapezoid, not triangle) condition, holds iff LHS > RHS.
function [LHS, RHS, isTrap] = nondegeneracyEq45(d, g, theta, amax)
    LHS = (d.Xf - d.Xr) * g * cos(theta);
    RHS = 2 * d.Zcm * amax;
    isTrap = LHS > RHS;
end

% Eq. (3.46)/(3.47), critical acceleration, i.e. Eq.(3.45) restated as a bound on a_max.
function acrit = criticalAccelEq46(d, g, theta)
    acrit = (d.Xf - d.Xr) * g * cos(theta) / (2 * d.Zcm);
end

% Eq. (3.50) antiderivative, W(Fz), the safe-region width, integrated.
% Psi over [a,b] (Eq. 3.51) = psiAntiderivative(...,b) - psiAntiderivative(...,a).
function P = psiAntiderivative(d, g, theta, amax, F)
    P = (d.Xf - d.Xr)/d.H * ( d.mr*g*cos(theta)*F + 0.5*F.^2 ) ...
        - 2*d.mr*d.Zcm*amax .* F / d.H;
end


%% Pipeline entry point
% This is what stability_metric.m calls (Plot=false), and what the manual
% input section calls too. Keep the accepted option names and the
% comOutsideWheelbase warning ID stable: stability_grid_pipeline.m silences
% that warning ID by name during parameter sweeps.

function result = coreRun(design, opts, name)

    validateGeometry(design, name);

    g     = opts.g;
    theta = deg2rad(opts.theta_deg);
    amax  = opts.a_max;
    d     = design;

    fn.Fz_back = @(Fx, ax) fzBack(Fx, ax, d, g, theta);
    fn.Fz_forw = @(Fx, ax) fzForw(Fx, ax, d, g, theta);
    fn.Fx_back = @(Fz, ax) fxBack(Fz, ax, d, g, theta);
    fn.Fx_forw = @(Fz, ax) fxForw(Fz, ax, d, g, theta);

    Fz_c                            = crossingPointFzc(d, g, theta, amax);      % Eq. (3.43)
    [~, ~, isTrap]                  = nondegeneracyEq45(d, g, theta, amax);     % Eq. (3.45)
    a_crit                          = criticalAccelEq46(d, g, theta);           % Eq. (3.46)/(3.47)

    Fzmax = opts.Fzmax;
    Fz0   = opts.FzFillMin;
    Psi        = psiAntiderivative(d,g,theta,amax,Fzmax) - psiAntiderivative(d,g,theta,amax,Fz0);  % Eq. (3.51)
    Psi_static = psiAntiderivative(d,g,theta,0,   Fzmax) - psiAntiderivative(d,g,theta,0,   Fz0);

    result = struct();
    result.name        = char(name);
    result.Fz_c        = Fz_c;
    result.isTrapezoid = isTrap;
    result.a_crit       = a_crit;
    result.Psi          = Psi;
    result.Psi_static   = Psi_static;
    result.fn            = fn;
    result.axes          = [];

    if opts.Verbose
        fprintf('%-14s  Fz,c=%+8.2f N (%s, Eq.3.43/3.45)   a_crit=%7.3f m/s^2 (Eq.3.46/47)   Psi=%9.1f N^2 (Eq.3.51)\n', ...
            name, Fz_c, tern(isTrap,'trapezoid','triangle '), a_crit, Psi);
    end

    if opts.Plot
        Fx_at_Fzc = fxBack(Fz_c, amax, d, g, theta);   % only needed to place the crossing-point marker
        result.axes = drawEnvelope(d, g, theta, amax, Fzmax, Fz0, opts.Visible, name, Fz_c, Fx_at_Fzc, isTrap, fn);
    end
end


function validateGeometry(d, name)
    req = ["mr","H","Xr","Xf","Xcm","Zcm"];
    for f = req
        if ~isfield(d, f) || ~isscalar(d.(f)) || ~isfinite(d.(f))
            error('stability_envelope_thesis:badGeometry', ...
                '%s: design.%s must be present and a finite scalar.', name, f);
        end
    end
    if ~(d.H > 0 && d.Xr < 0 && d.Xf > 0)
        error('stability_envelope_thesis:badGeometry', ...
            '%s: need H > 0, Xr < 0, Xf > 0 (got H=%.4f, Xr=%.4f, Xf=%.4f).', name, d.H, d.Xr, d.Xf);
    end
    if d.Xcm <= d.Xr || d.Xcm >= d.Xf
        warning('stability_envelope_thesis:comOutsideWheelbase', ...
            '%s: CoM Xcm=%.4f is outside the wheelbase [%.4f, %.4f].', name, d.Xcm, d.Xr, d.Xf);
    end
end


%% Manual input
% Not called by the pipeline. Edit the design struct and knobs below
% directly, then press Run (or call this file with no arguments) to plot
% that envelope. Opens a MATLAB figure, nothing written to disk.

function manualInputSection()

    % EDIT ME: illustrative geometry only, not a real concept
    design = struct( ...
        'mr',  19.624, ...     % total rollator mass [kg]
        'H',   1.075, ...      % handle height above O [m]
        'Xr', -0.300, ...      % rear wheel axis position rel. to O [m]  (Xr < 0)
        'Xf',  0.300, ...      % front wheel axis position rel. to O [m] (Xf > 0)
        'Xcm', 0.001672, ...   % CoM fore/aft position rel. to O [m]
        'Zcm', 0.211226);      % CoM height rel. to O [m]

    name      = "Example design";
    theta_deg = 0;      % ground incline [deg], front-downhill positive
    a_max     = 1.5;    % prescribed |a_x| bound [m/s^2]
    Fzmax     = 300;    % vertical force ceiling [N] (fill top + Psi upper limit)

    stability_envelope_thesis(design, 'Name', name, ...
        'theta_deg', theta_deg, 'a_max', a_max, 'Fzmax', Fzmax);

end


%% Plotting
% Presentation only, not thesis content. Styling fixed below rather than
% exposed as options since nothing in the pipeline varies it.

function ax = drawEnvelope(d, g, theta, amax, Fzmax, Fz0, visible, name, Fz_c, Fx_at_Fzc, isTrap, fn)

    COL_BACK   = [0 0.4470 0.7410];        % backward tipping
    COL_FORW   = [0.8500 0.3250 0.0980];   % forward tipping
    COL_FILL   = [0.60 0.88 0.60];
    FILL_ALPHA = 0.35;
    SOLID_LW   = 2;
    DASH_LW    = 1.5;
    FONT_NAME  = 'Times New Roman';
    FONT_SIZE  = 11;
    GRID_ALPHA = 0.15;

    span = max(Fzmax - Fz0, eps);
    yl = [Fz0 - 0.15*span, Fzmax + 0.15*span];
    Fz = linspace(yl(1), yl(2), 400);

    Fx_r0  = fn.Fx_back(Fz,  0);       % solid, a_x = 0
    Fx_f0  = fn.Fx_forw(Fz,  0);       % solid, a_x = 0
    Fx_rHi = fn.Fx_back(Fz,  amax);    % dashed, backward binds at a_x = +a_max  (Eq. 3.41)
    Fx_fLo = fn.Fx_forw(Fz, -amax);    % dashed, forward  binds at a_x = -a_max  (Eq. 3.42)

    fig = newLightFigure(sprintf('%s - Envelope', name), visible);
    ax  = axes(fig);
    hold(ax, 'on');

    % Green fill: region safe for every a_x in [-a_max, +a_max], clipped to
    % [FzFillMin, Fzmax]. Between Fx_fLo (left) and Fx_rHi (right) wherever
    % left <= right, which handles the triangle case cleanly: below Fz_c the
    % bounds have crossed and nothing gets filled.
    Fzf   = linspace(Fz0, Fzmax, 300);
    left  = fn.Fx_forw(Fzf, -amax);
    right = fn.Fx_back(Fzf,  amax);
    ok = left <= right;
    if any(ok)
        xs = [left(ok), fliplr(right(ok))];
        ys = [Fzf(ok),  fliplr(Fzf(ok))];
        fill(ax, xs, ys, COL_FILL, 'FaceAlpha', FILL_ALPHA, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    plot(ax, Fx_r0,  Fz, '-',  'Color', COL_BACK, 'LineWidth', SOLID_LW);
    plot(ax, Fx_f0,  Fz, '-',  'Color', COL_FORW, 'LineWidth', SOLID_LW);
    plot(ax, Fx_rHi, Fz, '--', 'Color', COL_BACK, 'LineWidth', DASH_LW);
    plot(ax, Fx_fLo, Fz, '--', 'Color', COL_FORW, 'LineWidth', DASH_LW);

    if Fz_c >= yl(1) && Fz_c <= yl(2)
        plot(ax, Fx_at_Fzc, Fz_c, 'k.', 'MarkerSize', 12, 'HandleVisibility', 'off');
        tag = tern(isTrap, '  F_{z,c} (< 0: trapezoid)', '  F_{z,c} (> 0: triangle apex)');
        text(ax, Fx_at_Fzc, Fz_c, tag, 'Color', 'k', 'FontName', FONT_NAME, ...
            'FontSize', FONT_SIZE-2, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    end

    allFx = [Fx_r0, Fx_f0, Fx_rHi, Fx_fLo];
    lo = min(allFx); hi = max(allFx);
    m  = 0.08 * (hi - lo);
    xlim(ax, [lo - m, hi + m]);
    ylim(ax, yl);

    grid(ax, 'on'); box(ax, 'on');
    ax.GridAlpha = GRID_ALPHA;
    ax.FontName  = FONT_NAME;
    ax.FontSize  = FONT_SIZE;
    ax.Color     = 'white';
    ax.XColor    = [0 0 0];
    ax.YColor    = [0 0 0];

    xlabel(ax, 'Horizontal handle force F_x (N), positive toward user', 'Color', 'k');
    ylabel(ax, 'Vertical handle force F_z (N), positive downward', 'Color', 'k');
    title(ax, sprintf('%s: Stability Envelope', name), 'FontName', FONT_NAME, 'Color', 'k');

    % Inline rotated labels instead of a legend box. Solid and dashed on the
    % same side sit only a few N apart (the a_x shift), so each label gets
    % its own anchor fraction plus a horizontal offset off its line: the
    % tipping labels sit just outside the solid boundary, the a_x labels
    % just inside the dashed one.
    addRotatedLineLabel(ax, Fx_r0,  Fz, 'backward tipping', COL_BACK, 0.72,  22);
    addRotatedLineLabel(ax, Fx_f0,  Fz, 'forward tipping',  COL_FORW, 0.72, -22);
    addRotatedLineLabel(ax, Fx_rHi, Fz, sprintf('a_x = +%.2g m/s2', amax), COL_BACK, 0.45, -22);
    addRotatedLineLabel(ax, Fx_fLo, Fz, sprintf('a_x = -%.2g m/s2', amax), COL_FORW, 0.45,  22);

    drawnow;
end


function fig = newLightFigure(figName, visible)
%NEWLIGHTFIGURE  Figure forced to the classic light theme. R2025a+ can
% otherwise inherit the desktop dark theme's pale grey axis/text colors
% even with figure Color forced white, isprop-guarded for older releases.
    if visible, vis = 'on'; else, vis = 'off'; end
    fig = figure('Name', figName, 'Visible', vis, 'Color', 'white');
    if isprop(fig, 'Theme')
        try
            fig.Theme = 'light';
        catch
        end
    end
end


function addRotatedLineLabel(ax, xdata, ydata, txt, color, anchorFrac, xOffPts)
%ADDROTATEDLINELABEL  Text label rotated to a line's on-screen slope, used
% in place of a legend entry. Uses the axes' actual points-per-data-unit
% scale in each direction since this figure doesn't use axis equal, so a
% naive atan2d on raw data would give the wrong angle. anchorFrac (0-1)
% sets where along the line the label sits; xOffPts shifts the anchor
% horizontally by that many points (signed) so labels on near-coincident
% lines clear each other.

    if nargin < 7, xOffPts = 0; end

    n = numel(xdata);
    i1 = max(1, round(anchorFrac*n));
    i2 = min(n, i1 + 5);
    if i2 == i1
        return
    end

    prevUnits = ax.Units;
    ax.Units = 'points';
    posPts = ax.Position;
    ax.Units = prevUnits;

    xl = xlim(ax); yl = ylim(ax);
    pxPerX = posPts(3) / diff(xl);
    pxPerY = posPts(4) / diff(yl);

    dx = (xdata(i2) - xdata(i1)) * pxPerX;
    dy = (ydata(i2) - ydata(i1)) * pxPerY;
    angleDeg = atan2d(dy, dx);

    xAnchor = xdata(i1) + xOffPts / pxPerX;
    yAnchor = ydata(i1);

    % normalise into (-90, 90], MATLAB doesn't auto-flip rotated text to stay upright
    if angleDeg > 90
        angleDeg = angleDeg - 180;
    elseif angleDeg <= -90
        angleDeg = angleDeg + 180;
    end

    % interpreter 'none': the default 'tex' turns "a_x" into a subscript
    % that visually detaches once the text is rotated to a steep angle
    text(ax, xAnchor, yAnchor, char(txt), 'Color', color, ...
        'Rotation', angleDeg, 'FontSize', 9, 'FontName', ax.FontName, ...
        'Interpreter', 'none', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end


function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end

function T = stability_grid_pipeline(opts)
%STABILITY_GRID_PIPELINE  (Xr, Xcm) grid search at a fixed wheelbase L, driven
% through stability_metric.m, with Excel in/out and a (Xr, Xcm) balance figure.
% Edit the arguments block below, then press Run.
%
% adoptable = feasible (Xr < Xcm < Xf) & heel_clearance_ok (Xr <= XrHeelLimit).
% Xr_hi is allowed to run a bit past XrHeelLimit on purpose, those rows are
% real evaluations, just not adoptable, so the balance figure can show what
% happens if the heel-clearance requirement moves.
%
% Heads up: min |Psi_r - Psi_f| does NOT have a unique minimiser, its zero
% set is a whole line (Xcm_balance below), and Psi itself is constant across
% the grid for fixed L. "Best" grid point = closest lattice node to that line.

arguments
    opts.L               (1,1) double  = 1.20      % fixed wheelbase (m)
    opts.mr              (1,1) double  = 32.964    % rollator mass (kg), CAD
    opts.H               (1,1) double  = 1.075     % handle height (m)
    opts.Zcm             (1,1) double  = 0.259840  % CoM height (m), CAD

    opts.Xr_lo           (1,1) double  = -1.19
    opts.Xr_hi           (1,1) double  = -0.70
    opts.XrHeelLimit     (1,1) double  = -0.85     % eq. 4.2 design requirement

    opts.DXr             (1,1) double  = 0.01
    opts.DXcm            (1,1) double  = 0.01

    opts.theta_deg       (1,1) double  = 0
    opts.a_max           (1,1) double  = 1.5
    opts.Fzmax           (1,1) double  = 300
    opts.FzFillMin       (1,1) double  = 0
    opts.g               (1,1) double  = 9.81
    opts.r               (1,1) double  = 2

    opts.Plot            (1,1) logical = true
    opts.OutDir          (1,1) string  = ""        % "" -> ./grid_results next to this file
    opts.InName          (1,1) string  = "stability_grid_in.xlsx"
    opts.OutName         (1,1) string  = "stability_grid_out.xlsx"
    opts.SavePlot        (1,1) string  = ""
    opts.Verbose         (1,1) logical = true
end

here = fileparts(mfilename('fullpath'));

% This file lives in "5 Optimization"; stability_metric.m and
% walking_force_ellipses.m live in the sibling folder "4 Stability Metric"
% (which in turn reaches "2 Stability Envelope" on its own), make sure
% they're reachable regardless of the caller's own folder / MATLAB path state.
if exist('stability_metric', 'file') ~= 2
    addpath(fullfile(here, '..', '4 Stability Metric'));
end

outDir = opts.OutDir;
if strlength(outDir) == 0
    % one self-contained folder per run: date + the two parameters that
    % define the design point, so grid_results/ accumulates a legible
    % history instead of overwriting the same two files every time.
    runName = string(datetime('now', 'Format', 'yyyy-MM-dd')) + ...
        sprintf('_L%.2f_XrHeel%.2f', opts.L, opts.XrHeelLimit);
    outDir = string(fullfile(here, 'grid_results', runName));
end
if ~isfolder(outDir); mkdir(outDir); end
inPath  = char(fullfile(outDir, opts.InName));
outPath = char(fullfile(outDir, opts.OutName));

Xr_vals = round((opts.Xr_lo : opts.DXr : opts.Xr_hi).', 10);
d_axis  = round((0 : opts.DXcm : opts.L).', 10);
[XR, D] = ndgrid(Xr_vals, d_axis);
XCM     = XR + D;
XR = XR(:);  XCM = XCM(:);
Xf = XR + opts.L;
n  = numel(XR);

feasible       = (XCM > XR) & (XCM < Xf);
heel_clearance_ok = XR <= opts.XrHeelLimit + 1e-12;
adoptable      = feasible & heel_clearance_ok;

Gin = table( ...
    (1:n).', ...
    repmat(opts.mr,n,1), repmat(opts.H,n,1), XR, Xf, XCM, repmat(opts.Zcm,n,1), ...
    repmat(opts.L,n,1), ...
    repmat(opts.theta_deg,n,1), repmat(opts.a_max,n,1), repmat(opts.Fzmax,n,1), ...
    repmat(opts.FzFillMin,n,1), repmat(opts.g,n,1), repmat(opts.r,n,1), ...
    repmat(opts.XrHeelLimit,n,1), ...
    feasible, heel_clearance_ok, adoptable, ...
    'VariableNames', {'designID','mr','H','Xr','Xf','Xcm','Zcm', 'L', ...
        'theta_deg','a_max','Fzmax','FzFillMin','g', ...
        'r_sigma','XrHeelLimit', ...
        'feasible','heel_clearance_ok','adoptable'});

if isfile(inPath); delete(inPath); end
writetable(Gin, inPath);
if opts.Verbose
    fprintf('Grid: %d Xr x %d Xcm = %d designs   (%d feasible, %d adoptable)  ->  %s\n', ...
        numel(Xr_vals), n/numel(Xr_vals), n, sum(feasible), sum(adoptable), inPath);
end

G     = readtable(inPath, 'TextType', 'string');
nrows = height(G);
feasible  = logical(G.feasible);
adoptable = logical(G.adoptable);

E          = walking_force_ellipses('Verbose', false);
K          = numel(E);
condLabels = strings(1,K);
for k = 1:K; condLabels(k) = E(k).label; end

z = nan(nrows,1);
Psi=z; Psi_r=z; Psi_f=z; Psi_diff=z; absPsi_diff=z; Psi_ratio=z;
Psi_r_closed=z; Psi_f_closed=z; R_margin=z; F_margin=z; Fz_c=z; a_crit=z;
clampFree    = false(nrows,1);
validR       = false(nrows,1);
validF       = false(nrows,1);
isTrapezoid  = false(nrows,1);
com_on_axle  = false(nrows,1);
d_mu   = nan(nrows,K);   d_ell = nan(nrows,K);
ell_bind = strings(nrows,K);
crossAny = false(nrows,K);
d_mu_min=z; d_ell_min=z;

ws          = warning('off', 'stability_envelope_thesis:comOutsideWheelbase');
restoreWarn = onCleanup(@() warning(ws));

t0 = tic;
for i = 1:nrows
    d = struct('mr',G.mr(i), 'H',G.H(i), 'Xr',G.Xr(i), 'Xf',G.Xf(i), ...
               'Xcm',G.Xcm(i), 'Zcm',G.Zcm(i));
    com_on_axle(i) = (d.Xcm <= d.Xr) || (d.Xcm >= d.Xf);

    m = stability_metric(d, E, 'Name', "g"+i, ...
        'theta_deg', G.theta_deg(i), 'a_max', G.a_max(i), 'g', G.g(i), ...
        'Fzmax', G.Fzmax(i), 'FzFillMin', G.FzFillMin(i), 'r', G.r_sigma(i), ...
        'Verbose', false);

    Psi(i)=m.Psi;  Psi_r(i)=m.Psi_r;  Psi_f(i)=m.Psi_f;
    Psi_diff(i)=m.Psi_diff;  absPsi_diff(i)=abs(m.Psi_diff);  Psi_ratio(i)=m.Psi_ratio;
    Psi_r_closed(i)=m.split.Psi_r_closed;  Psi_f_closed(i)=m.split.Psi_f_closed;
    R_margin(i)=m.split.R_LHS - m.split.R_RHS;
    F_margin(i)=m.split.F_LHS - m.split.F_RHS;
    clampFree(i)=m.split.clampFree;  validR(i)=m.split.validR;  validF(i)=m.split.validF;
    isTrapezoid(i)=m.isTrapezoid;  Fz_c(i)=m.Fz_c;  a_crit(i)=m.a_crit;
    for k = 1:K
        d_mu(i,k)=m.cond(k).d_mu;   d_ell(i,k)=m.cond(k).d_ellipse;
        ell_bind(i,k)=m.cond(k).ellipse_binding;
        crossAny(i,k)=m.cond(k).crosses_backward || m.cond(k).crosses_forward;
    end
    d_mu_min(i)=m.d_mu_min;  d_ell_min(i)=m.d_ellipse_min;

    if opts.Verbose && mod(i,500)==0
        fprintf('  %5d / %d   (%.0f s)\n', i, nrows, toc(t0));
    end
end

R = G;
R.Psi            = Psi;
R.Psi_r          = Psi_r;
R.Psi_f          = Psi_f;
R.Psi_diff       = Psi_diff;
R.absPsi_diff    = absPsi_diff;
R.Psi_ratio      = Psi_ratio;
R.Psi_r_closed   = Psi_r_closed;
R.Psi_f_closed   = Psi_f_closed;
R.clampFree      = clampFree;
R.validR         = validR;
R.validF         = validF;
R.R_margin       = R_margin;
R.F_margin       = F_margin;
R.isTrapezoid_Eq20 = isTrapezoid;
R.Fz_c           = Fz_c;
R.a_crit         = a_crit;
R.com_on_axle    = com_on_axle;
for k = 1:K
    R.("d_mu_c"  + k) = d_mu(:,k);
    R.("d_ell_c" + k) = d_ell(:,k);
    R.("bind_c"  + k) = ell_bind(:,k);
    R.("cross_c" + k) = crossAny(:,k);
end
R.d_mu_min  = d_mu_min;
R.d_ell_min = d_ell_min;

% analytic exact-balance line  Psi_r = Psi_f  <=>  Xcm = Xcm_balance
Xcm_balance           = (G.Xr + G.L/2) .* (1 + G.Fzmax ./ (2 .* G.mr .* G.g));
R.Xcm_balance         = Xcm_balance;
R.offset_from_balance = G.Xcm - Xcm_balance;

% rank by the objective over ADOPTABLE rows only (1 = smallest |Psi_r - Psi_f|)
absF = absPsi_diff;  absF(~adoptable) = NaN;
[~, ord] = sort(absF, 'ascend', 'MissingPlacement', 'last');
rankv = nan(nrows,1);
na = sum(adoptable & ~isnan(absPsi_diff));
rankv(ord(1:na)) = (1:na).';
R.rank_absPsi_diff_adoptable = rankv;

if isfile(outPath); delete(outPath); end
writetable(R, outPath);

% sanity checks + console summary
psiSpread = max(Psi) - min(Psi);
psiConst  = psiSpread <= 1e-6 * max(abs(mean(Psi)), 1);
nTrapFail = sum(~isTrapezoid);
nNotClamp = sum(~clampFree(adoptable));
clampCM   = 100 * mean(G.Zcm) * mean(G.a_max) / mean(G.g);
[bestAbs, iBest] = min(absF);
onLine = adoptable & abs(R.offset_from_balance) <= opts.DXcm/2 + 1e-9;

if opts.Verbose
    fprintf('\n==================== grid summary ====================\n');
    fprintf('  designs run              : %d   (%d feasible, %d adoptable)\n', ...
        nrows, sum(feasible), sum(adoptable));
    fprintf('  Psi (must be constant)   : %.4f N^2   spread %.2e N^2   %s\n', ...
        mean(Psi), psiSpread, tern(psiConst, 'constant, ok', ...
        '*** NOT CONSTANT, INVESTIGATE ***'));
    fprintf('  Eq.20 trapezoid holds    : %d / %d rows\n', nrows - nTrapFail, nrows);
    fprintf('  Psi split clamp-free     : %d / %d adoptable rows  (rest: CoM within ~%.1f cm of an axle)\n', ...
        sum(adoptable) - nNotClamp, sum(adoptable), clampCM);
    fprintf('  objective: min |Psi_r - Psi_f| over adoptable rows\n');
    fprintf('  best adoptable grid point : Xr* = %+.4f m   Xcm* = %+.4f m\n', ...
        G.Xr(iBest), G.Xcm(iBest));
    fprintf('                            : |Psi_r - Psi_f| = %.2f N^2   (Psi_r %.1f, Psi_f %.1f, Psi %.1f)\n', ...
        bestAbs, Psi_r(iBest), Psi_f(iBest), Psi(iBest));
    if any(onLine)
        fprintf('  adoptable points ON the balance line (|offset| <= %.0f mm): %d\n', ...
            1000*opts.DXcm/2, sum(onLine));
        fprintf('     Xr  in [%+.3f, %+.3f] m   Xcm in [%+.3f, %+.3f] m\n', ...
            min(R.Xr(onLine)), max(R.Xr(onLine)), min(R.Xcm(onLine)), max(R.Xcm(onLine)));
    end
    fprintf('  balance line : Xcm*(Xr) = (Xr + L/2)(1 + Fzmax/(2 mr g)) = %.4f*Xr + %.4f\n', ...
        1 + opts.Fzmax/(2*opts.mr*opts.g), (opts.L/2)*(1 + opts.Fzmax/(2*opts.mr*opts.g)));
    fprintf('  ellipse conditions: c1 = %s\n', condLabels(1));
    for k = 2:K
        fprintf('                     c%d = %s\n', k, condLabels(k));
    end
    fprintf('  workbooks:\n    %s\n    %s\n', inPath, outPath);
    fprintf('=====================================================\n');
end

if ~psiConst
    warning('stability_grid_pipeline:psiNotConstant', ...
        'Psi varies by %.3e N^2 across the grid, it should be constant for fixed L.', psiSpread);
end

% 5. balance figure
if opts.Plot
    plotBase = opts.SavePlot;
    if strlength(plotBase) == 0
        plotBase = string(fullfile(outDir, "stability_grid_balance_map"));
    end
    savePaths = [plotBase + ".svg", plotBase + ".png"];  % TODO: add .pdf here too if the appendix ends up wanting it
    pad = max(opts.DXr, opts.DXcm);
    xrLim = [min(R.Xr) - pad, max(R.Xr) + pad];
    % centred on the balance line, not the raw Xcm sweep extent (too wide, unreadable aspect ratio)
    balMargin = 0.12;
    kbal = 1 + opts.Fzmax / (2 * opts.mr * opts.g);
    xcmBalEnds = kbal * (xrLim + opts.L/2);
    xcmTop = -0.20;
    xcmLim = [min(xcmBalEnds) - balMargin, xcmTop];
    try
        stability_balance_map(R, 'Save', savePaths, 'Xr_lim', xrLim, 'Xcm_lim', xcmLim, ...
            'R26', opts.XrHeelLimit);
        if opts.Verbose
            fprintf('  figure: %s  (+ .png)\n', savePaths(1));
        end
    catch ME
        warning('stability_grid_pipeline:plotFailed', ...
            'Balance figure failed: %s', ME.message);
    end
end

if nargout > 0
    T = R;
end

end % stability_grid_pipeline


function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end

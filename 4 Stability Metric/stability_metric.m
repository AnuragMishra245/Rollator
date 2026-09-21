function m = stability_metric(design, E, opts)
%STABILITY_METRIC  Psi (Eq. 3.51), its rearward/forward split (Sec.
% 3.4.3/3.4.4), and ellipse-to-boundary distances (Sec. 3.4.5) for one design.
%
% Press Run (or call with no args) for the manual input section. Called
% from the pipeline with a design struct, normal path below, silent unless
% Verbose.
%
%   E = walking_force_ellipses();      % once
%   for i = 1:numel(designs)
%       m(i) = stability_metric(designs(i), E);
%   end
%
% design needs mr, H, Xr, Xf, Xcm, Zcm (Xr < 0 < Xf). E is optional, defaults
% to walking_force_ellipses() with its default folder. opts.a_max = 0
% measures against the static boundary instead.
%
% Sign convention: d > 0 = safe side of that edge, d < 0 for the ellipse
% means it pokes through by |d|.

arguments
    design = []                     % [] => manual input section (edit + press Run)
    E = []
    opts.Name       (1,1) string  = "design"
    opts.r          (1,1) double  = 2
    opts.theta_deg  (1,1) double  = 0
    opts.a_max      (1,1) double  = 1.5
    opts.g          (1,1) double  = 9.81
    opts.Fzmax      (1,1) double  = 300
    opts.FzFillMin  (1,1) double  = 0
    opts.Verbose    (1,1) logical = true
end

% This file lives in "4 Stability Metric"; stability_envelope_thesis.m lives
% in the sibling folder "2 Stability Envelope", make sure it's reachable
% regardless of the caller's own folder / MATLAB path state.
if exist('stability_envelope_thesis', 'file') ~= 2
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', '2 Stability Envelope'));
end

if isempty(design) || ~isstruct(design) || isempty(fieldnames(design))
    manualInputSection();
    m = [];
    return
end
validateDesign(design);

if isempty(E)
    E = walking_force_ellipses('Verbose', false);
end

g   = opts.g;
th  = deg2rad(opts.theta_deg);
r   = opts.r;
amx = opts.a_max;
d   = design;

% Psi and envelope scalars (no plot, no asserts, fast in a sweep)
env = stability_envelope_thesis(d, 'Name', opts.Name, ...
        'theta_deg', opts.theta_deg, 'a_max', amx, 'g', g, ...
        'Fzmax', opts.Fzmax, 'FzFillMin', opts.FzFillMin, ...
        'Plot', false, 'Verbose', false, 'RunAsserts', false);

% Eq. (3.28)/(3.29) constants, shared by the margin lines below and by
% psiSplit's edge functions, computed once here rather than re-derived
% in each place that needs them.
Kb = d.mr*g*cos(th)*(d.Xcm - d.Xr) + d.mr*g*d.Zcm*sin(th);
Kf = d.mr*g*cos(th)*(d.Xcm - d.Xf) + d.mr*g*d.Zcm*sin(th);

% Step 1b: rearward / forward split of Psi at Fx = 0
split = psiSplit(d, Kb, Kf, amx, opts.Fzmax, opts.FzFillMin);
if split.clampFree
    assert(abs((split.Psi_r + split.Psi_f) - env.Psi) < 1e-6*max(abs(env.Psi), 1), ...
        'stability_metric:psiSplitSum', ...
        'Psi_r + Psi_f = %.4f != Psi = %.4f (Eq. 3.51 identity broken).', ...
        split.Psi_r + split.Psi_f, env.Psi);
end

% margin-boundary lines  a_i*Fx + Fz - c_i = 0  (Eqs. 3.59/3.60)
% From Eqs. (3.28)/(3.29):  Fz = ( K_i - H*Fx - mr*Zcm*ax ) / X_i.
% Rearranged:  a_i = H / X_i ,  c_i = ( K_i - mr*Zcm*ax_edge ) / X_i .
% Backward edge binds at ax = +a_max, forward edge at ax = -a_max (the same
% edges stability_envelope_thesis fills the safe region between).
a_back = d.H / d.Xr;
c_back = (Kb - d.mr*d.Zcm*(+amx)) / d.Xr;
a_forw = d.H / d.Xf;
c_forw = (Kf - d.mr*d.Zcm*(-amx)) / d.Xf;

nrm_back = hypot(a_back, 1);           % |[a_i, 1]|
nrm_forw = hypot(a_forw, 1);
nhat_back = [a_back; 1] / nrm_back;    % unit outward normal of each edge
nhat_forw = [a_forw; 1] / nrm_forw;

% per-condition distances
K = numel(E);
condCell = cell(1, K);
tblRows = cell(K, 18);

for k = 1:K
    mu = E(k).mu(:);
    S  = E(k).S;

    d_back_mu = (a_back*mu(1) + mu(2) - c_back) / nrm_back;   % Eq. 3.61
    d_forw_mu = (a_forw*mu(1) + mu(2) - c_forw) / nrm_forw;

    w_back = r * sqrt(nhat_back.' * S * nhat_back);           % Eq. 3.63
    w_forw = r * sqrt(nhat_forw.' * S * nhat_forw);

    d_back_ell = d_back_mu - w_back;                          % Eq. 3.64
    d_forw_ell = d_forw_mu - w_forw;

    [d_mu,  iMu ] = min([d_back_mu,  d_forw_mu]);             % Eq. 3.62
    [d_ell, iEll] = min([d_back_ell, d_forw_ell]);            % Eq. 3.65
    edges = ["backward", "forward"];

    c = struct();
    c.label            = E(k).label;
    c.d_back_mu        = d_back_mu;    % [N]
    c.d_forw_mu        = d_forw_mu;
    c.d_mu             = d_mu;
    c.mu_binding       = edges(iMu);
    c.w_back           = w_back;
    c.w_forw           = w_forw;
    c.d_back_ell       = d_back_ell;
    c.d_forw_ell       = d_forw_ell;
    c.d_ellipse        = d_ell;
    c.ellipse_binding  = edges(iEll);
    c.crosses_backward = d_back_mu <= w_back;
    c.crosses_forward  = d_forw_mu <= w_forw;
    condCell{k} = c;

    tblRows(k,:) = {string(opts.Name), string(E(k).label), env.Psi, ...
        split.Psi_r, split.Psi_f, split.Psi_diff, ...
        d_mu, d_ell, d_back_mu, d_forw_mu, d_back_ell, d_forw_ell, ...
        w_back, w_forw, c.mu_binding, c.ellipse_binding, ...
        c.crosses_backward, c.crosses_forward};
end

cond = [condCell{:}];

m = struct();
m.name        = char(opts.Name);
m.Psi         = env.Psi;
m.Psi_static  = env.Psi_static;           % at a_x = 0
m.Psi_r       = split.Psi_r;              % rearward (pull) share, Fx > 0
m.Psi_f       = split.Psi_f;              % forward (push) share, Fx < 0
m.Psi_diff    = split.Psi_diff;
m.Psi_ratio   = split.Psi_ratio;          % 0.5 = balanced
m.split = rmfield(split, {'Psi_r','Psi_f','Psi_diff','Psi_ratio'});
m.Fz_c        = env.Fz_c;
m.isTrapezoid = env.isTrapezoid;
m.a_crit      = env.a_crit;
m.cond        = cond;
m.d_mu_min      = min([cond.d_mu]);
m.d_ellipse_min = min([cond.d_ellipse]);

m.tbl = cell2table(tblRows, 'VariableNames', ...
    {'Concept','Condition','Psi','Psi_r','Psi_f','Psi_diff', ...
     'd_mu','d_ellipse','d_back_mu','d_forw_mu', ...
     'd_back_ell','d_forw_ell','w_back','w_forw','MuBinding','EllipseBinding', ...
     'CrossesBack','CrossesForw'});         % one row per condition; Psi* repeated per design

if opts.Verbose
    report(m, opts);
end

end % stability_metric


function manualInputSection()
%MANUALINPUTSECTION  Not called by the pipeline. Edit the design struct and
% knobs below directly, then press Run (or call stability_metric with no
% arguments) to compute and print the stability metrics for that geometry,
% prints to the Command Window, nothing is written to disk.

    % EDIT ME: illustrative geometry only, not a real concept
    design = struct( ...
        'mr',  19.624, ...     % total rollator mass [kg]
        'H',   1.075, ...      % handle height above O [m]
        'Xr', -0.300, ...      % rear wheel axis position rel. to O [m]  (Xr < 0)
        'Xf',  0.300, ...      % front wheel axis position rel. to O [m] (Xf > 0)
        'Xcm', 0.001672, ...   % CoM fore/aft position rel. to O [m]
        'Zcm', 0.211226);      % CoM height rel. to O [m]

    name      = "Example design";
    r         = 2;      % ellipse scale [standard deviations]
    theta_deg = 0;       % ground incline [deg], front-downhill positive
    a_max     = 1.5;     % prescribed |a_x| bound [m/s^2]
    Fzmax     = 300;      % vertical force ceiling [N] (fill top + Psi upper limit)

    stability_metric(design, [], 'Name', name, 'r', r, ...
        'theta_deg', theta_deg, 'a_max', a_max, 'Fzmax', Fzmax);

end


function validateDesign(d)
    if ~isscalar(d)
        error('stability_metric:notScalar', 'Pass one design struct at a time.');
    end
    req = ["mr","H","Xr","Xf","Xcm","Zcm"];
    for f = req
        if ~isfield(d,f) || ~isscalar(d.(f)) || ~isfinite(d.(f))
            error('stability_metric:badDesign', 'design.%s must be a finite scalar.', f);
        end
    end
    if ~(d.H > 0 && d.Xr < 0 && d.Xf > 0)
        error('stability_metric:badGeometry', 'Need H > 0, Xr < 0, Xf > 0.');
    end
end


function s = psiSplit(d, Kb, Kf, amx, Fzmax, Fz0)
%PSISPLIT  Rearward/forward split of Psi about Fx = 0 (Sec. 3.4.3/3.4.4).
% Psi_r+Psi_f=Psi only when neither edge crosses Fx=0 over [Fz0,Fzmax] (the
% (R)/(F) conditions below) - numerical Psi_r/Psi_f are always valid, the
% closed forms are just for convenience and only match when (R)/(F) hold.

    mz = d.mr * d.Zcm;
    Fxr = @(Fz) (Kb - mz*amx - d.Xr.*Fz) / d.H;   % backward edge at a_x = +a_max
    Fxf = @(Fz) (Kf + mz*amx - d.Xf.*Fz) / d.H;   % forward  edge at a_x = -a_max

    % numerical (clamped): always valid. Waypoint at the sign change of
    %     each integrand, if it falls inside the range, so the adaptive
    %     quadrature resolves the kink exactly.
    brkR = (Kb - mz*amx) / d.Xr;                   % Fz where Fx_r = 0
    brkF = (Kf + mz*amx) / d.Xf;                   % Fz where Fx_f = 0
    wpR  = brkR(brkR > Fz0 & brkR < Fzmax);
    wpF  = brkF(brkF > Fz0 & brkF < Fzmax);
    Psi_r_num = integral(@(Fz) max( Fxr(Fz), 0), Fz0, Fzmax, 'Waypoints', wpR);
    Psi_f_num = integral(@(Fz) max(-Fxf(Fz), 0), Fz0, Fzmax, 'Waypoints', wpF);

    % closed form (clamp-free antiderivatives, generalised to Fz0)
    Ar = @(F) ( (Kb - mz*amx).*F - 0.5*d.Xr.*F.^2 ) / d.H;
    Af = @(F) ( -(Kf + mz*amx).*F + 0.5*d.Xf.*F.^2 ) / d.H;
    Psi_r_cf = Ar(Fzmax) - Ar(Fz0);
    Psi_f_cf = Af(Fzmax) - Af(Fz0);

    % clamp-free conditions (R)/(F): integrand >= 0 at both ends
    tol    = 1e-9 * (abs(Kb) + abs(Kf) + 1);
    validR = Fxr(Fz0) >= -tol && Fxr(Fzmax) >= -tol;
    validF = -Fxf(Fz0) >= -tol && -Fxf(Fzmax) >= -tol;

    s = struct();
    s.Psi_r        = Psi_r_num;
    s.Psi_f        = Psi_f_num;
    s.Psi_diff     = Psi_r_num - Psi_f_num;                        % Eq. 3.57
    s.Psi_ratio    = Psi_r_num / (Psi_r_num + Psi_f_num);          % Eq. 3.58
    s.Psi_r_closed = Psi_r_cf;
    s.Psi_f_closed = Psi_f_cf;
    s.clampFree    = validR && validF;
    s.validR       = validR;
    s.validF       = validF;
    s.R_LHS = Kb;               % Eq. 3.55 LHS (Fz0 = 0 form)
    s.R_RHS = mz*amx;
    s.F_LHS = -Kf;              % Eq. 3.56 LHS
    s.F_RHS = mz*amx;
end


function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end


function report(m, opts)
    fprintf('\n=== %s : stability metric ===\n', m.name);
    fprintf('  Psi = %.1f N^2   Psi_static = %.1f N^2   Fz,c = %+.2f N   trapezoid: %s\n', ...
        m.Psi, m.Psi_static, m.Fz_c, string(m.isTrapezoid));
    fprintf('  Psi_r = %.1f N^2 (pull, Fx>0)   Psi_f = %.1f N^2 (push, Fx<0)   Psi_r - Psi_f = %+.1f N^2   Psi_r/(Psi_r+Psi_f) = %.3f\n', ...
        m.Psi_r, m.Psi_f, m.Psi_diff, m.Psi_ratio);
    if m.split.clampFree
        fprintf('    (R)/(F) hold: closed form applies, Psi_r + Psi_f = Psi to %.1e N^2\n', ...
            abs(m.Psi_r + m.Psi_f - m.Psi));
    else
        fprintf(['    [!] Psi split NOT clamp-free: (R) %s (%.4g vs %.4g), (F) %s (%.4g vs %.4g).\n', ...
                 '        Psi_r/Psi_f above are from clamped numerical integration; the boxed\n', ...
                 '        closed form and the identity Psi_r + Psi_f = Psi do NOT apply here.\n'], ...
            tern(m.split.validR,'ok','VIOLATED'), m.split.R_LHS, m.split.R_RHS, ...
            tern(m.split.validF,'ok','VIOLATED'), m.split.F_LHS, m.split.F_RHS);
    end
    fprintf('  margin boundary: theta = %.3g deg,  a_x edges = +/- %.3g m/s^2,  ellipse r = %.3g sigma\n', ...
        opts.theta_deg, opts.a_max, opts.r);
    fprintf('  %-30s  %9s %9s  %9s %9s  %-9s %-9s\n', ...
        'condition', 'd(mu)', 'd(ell)', 'd_back_mu', 'd_forw_mu', 'muBind', 'ellBind');
    for k = 1:numel(m.cond)
        c = m.cond(k);
        flag = "";
        if c.crosses_backward || c.crosses_forward, flag = "  [ellipse crosses an edge]"; end
        fprintf('  %-30s  %9.2f %9.2f  %9.2f %9.2f  %-9s %-9s%s\n', ...
            c.label, c.d_mu, c.d_ellipse, c.d_back_mu, c.d_forw_mu, ...
            c.mu_binding, c.ellipse_binding, flag);
    end
    fprintf('  worst over conditions:  d(mu)_min = %.2f N   d(ellipse)_min = %.2f N\n', ...
        m.d_mu_min, m.d_ellipse_min);
end

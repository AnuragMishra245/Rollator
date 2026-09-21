function E = walking_force_ellipses(dataDir, opts)
%WALKING_FORCE_ELLIPSES  Per-condition (mu, S) of the measured walking handle
% forces, in the stability model's (Fx, Fz) convention. Computed once and
% reused for every candidate design in a stability_metric sweep.
%
% Returns a 1x3 struct array (.label, .n, .mu, .S).
%
% TREATMENT (frozen, identical to loadWalkingTrajectories.m and to
% interaction_forces_analysis.ipynb / covariance_ellipse.ipynb):
%   raw CSV columns  t_s, Fx_left, Fy_left, Fz_left, Fx_right, Fy_right, Fz_right
%   are in kilogram-force (kgf), sensor-local axes.
%       model Fz(t) =  9.81 * (Fz_left + Fz_right)      vertical, positive down
%       model Fx(t) = -9.81 * (Fy_left + Fy_right)      fore-aft, sign flipped
%   The sensor's own left/right "Fx" is mediolateral and isn't used here.
%   Nothing trimmed, filtered, or subsampled.
%
% Every load is checked against the handover reference means; a mismatch
% stops execution rather than returning silently wrong (Fx, Fz).

arguments
    dataDir (1,1) string = defaultDataDir()
    opts.Verbose (1,1) logical = true
end

g = 9.81;

cond = struct( ...
    'label', {'Constant velocity, v = 0.2 m/s', ...
              'Admittance, dv = 20, dw = 10', ...
              'Admittance, dv = 200, dw = 100'}, ...
    'stem',  {'walker_straight_1_constvel-0.2ms', ...
              'walker_straight_2_adm-dv20-dw10', ...
              'walker_straight_3_adm-dv200-dw100'}, ...
    'refFx', {-24.3, -31.4, -77.6}, ...     % handover reference means [N]
    'refFz', {262.3, 250.8, 253.2});

if ~isfolder(dataDir)
    error('walking_force_ellipses:missingDir', ...
        'Data folder not found: %s', dataDir);
end

E = repmat(struct('label', "", 'n', 0, 'mu', [], 'S', []), 1, numel(cond));

for k = 1:numel(cond)
    fpath = fullfile(dataDir, cond(k).stem + "_force.csv");
    if ~isfile(fpath)
        error('walking_force_ellipses:missingFile', 'Force file not found: %s', fpath);
    end

    F  = readtable(fpath);
    Fz =  g * (F.Fz_left + F.Fz_right);
    Fx = -g * (F.Fy_left + F.Fy_right);

    mu = [mean(Fx); mean(Fz)];
    S  = cov([Fx Fz]);          % 2x2, normalised by N-1

    dFx = abs(mu(1) - cond(k).refFx);
    dFz = abs(mu(2) - cond(k).refFz);
    assert(dFx < 1 && dFz < 1, 'walking_force_ellipses:refMismatch', ...
        ['%s: mean (Fx=%.2f, Fz=%.2f) does not match the reference table ', ...
         '(ref Fx=%.1f, Fz=%.1f; diff %.2f / %.2f), treatment or axis mapping is wrong.'], ...
        cond(k).label, mu(1), mu(2), cond(k).refFx, cond(k).refFz, dFx, dFz);

    E(k).label = string(cond(k).label);
    E(k).n     = height(F);
    E(k).mu    = mu;
    E(k).S     = S;

    if opts.Verbose
        fprintf(['  %-32s n=%6d | mu = [%+7.2f, %+7.2f] N | ', ...
                 'S = [%8.2f %9.2f ; %9.2f %9.2f] N^2\n'], ...
            cond(k).label, E(k).n, mu(1), mu(2), S(1,1), S(1,2), S(2,1), S(2,2));
    end
end

end % walking_force_ellipses


function d = defaultDataDir()
% CSVs live in "1 Interaction Forces Experiment/Interaction Force Measured
% Data", a sibling of this file's folder ("4 Stability Metric") under the
% "Stability Analysis" root (mirrors covariance_ellipse.ipynb's DATA_DIR).
    here = fileparts(mfilename('fullpath'));
    d = string(fullfile(here, '..', '1 Interaction Forces Experiment', 'Interaction Force Measured Data'));
end

clear
close all
clc

% Differentiates the filtered velocity (vx_filt, from velocity_conditioning.m)
% to get a(t) = gradient(vx_filt, t_s). a(t) itself is not filtered, this is
% just for looking at it. Phase segmentation (accel/steady/decel) isn't done
% here yet, want to look at the plots first.

conditions = {
    'const_v', 'Constant velocity 0.2 m/s';
    'adm1',    'Admittance 1 (dv=20, dw=10)';
    'adm2',    'Admittance 2 (dv=200, dw=100)'
    };

nCond = size(conditions, 1);
scriptDir = fileparts(mfilename('fullpath'));
resultsRoot = fullfile(scriptDir, 'Filtered Data');  % written by velocity_conditioning.m

edgeWindowSec = 8;   % width of start/end transient zoom windows, s
midWindowSec  = 8;   % width of the middle steady-region zoom window, s

peakCond = cell(nCond*2, 1);
peakSignal = cell(nCond*2, 1);
peakValue = nan(nCond*2, 1);

% Rough steady-state window, NOT the real accel/steady/decel segmentation.
% Just collecting a crude middle-of-trial mean/std of a(t) to justify
% treating a_x ~ 0 in the stability analysis for now.
prov_trimFracs   = [0.25, 0.20, 0.30];          % drop this fraction off each end
prov_midPct      = 100 * (1 - 2*prov_trimFracs); % -> 50%, 60%, 40% kept
prov_nWin        = numel(prov_trimFracs);
prov_rows        = nCond * prov_nWin;
prov_trial       = cell(prov_rows, 1);
prov_windowPct   = nan(prov_rows, 1);
prov_meanAx      = nan(prov_rows, 1);
prov_stdAx       = nan(prov_rows, 1);
prov_maxAbsAx    = nan(prov_rows, 1);
prov_r           = 0;

for k = 1:nCond
    condName  = conditions{k,1};
    condLabel = conditions{k,2};

    outDir = fullfile(resultsRoot, condName);
    srcFile = fullfile(outDir, 'filtered_velocity.csv');

    fprintf('\n=== %s (%s) ===\n', condLabel, condName);

    if ~exist(srcFile, 'file')
        error('acceleration_analysis:MissingInput', ...
            '%s: filtered_velocity.csv not found in %s', condName, outDir);
    end

    % load filtered velocity
    T = readtable(srcFile, 'Delimiter', ',');
    t_s = T.t_s(:);
    vx_filt = T.vx_filt(:);
    vx_raw = T.vx_raw(:);

    if any(isnan(t_s)) || any(isnan(vx_filt)) || any(isnan(vx_raw))
        error('acceleration_analysis:NaN', '%s: NaNs found in filtered_velocity.csv', condName);
    end
    if ~all(diff(t_s) > 0)
        error('acceleration_analysis:TimeNotIncreasing', '%s: t_s is not strictly increasing', condName);
    end

    % a_raw uses the same central-difference gradient but on the unfiltered
    % velocity - kept around so we can tell whether ringing near sharp
    % transitions in a(t) is a filtfilt artifact or was already in the raw signal.
    a = gradient(vx_filt, t_s);
    a_raw = gradient(vx_raw, t_s);

    outTable = table(t_s, vx_filt, a, a_raw, 'VariableNames', {'t_s', 'vx_filt', 'a', 'a_raw'});
    writetable(outTable, fullfile(outDir, 'acceleration_data.csv'))

    % descriptive summary, no phase detection yet
    a_min = min(a);
    a_max = max(a);
    a_mean = mean(a);
    fprintf('%s: min(a) = %.4f m/s^2, max(a) = %.4f m/s^2, mean(a) = %.4f m/s^2\n', ...
        condName, a_min, a_max, a_mean);

    % full-trial peak magnitude, filtered vs unfiltered - feeds peak_accel_table.csv below
    peak_a = max(abs(a));
    peak_a_raw = max(abs(a_raw));

    row = (k-1)*2 + 1;
    peakCond{row}     = condLabel;
    peakSignal{row}   = 'a (filtered)';
    peakValue(row)    = peak_a;
    peakCond{row+1}   = condLabel;
    peakSignal{row+1} = 'a_raw (unfiltered)';
    peakValue(row+1)  = peak_a_raw;

    fprintf('%s: full-trial peak |a| = %.4f m/s^2 (filtered), |a_raw| = %.4f m/s^2 (unfiltered)\n', ...
        condName, peak_a, peak_a_raw);

    t0 = t_s(1);
    t1 = t_s(end);
    duration = t1 - t0;

    % full-trial plot
    figFull = newfig();
    plot(t_s, a, 'Color', [0.85 0.325 0.098], 'LineWidth', 0.8)
    grid on
    box on
    xlabel('t (s)')
    ylabel('a (m/s^2)')
    title(sprintf('%s: Acceleration a(t), full trial', condLabel))
    saveas(figFull, fullfile(outDir, 'acceleration_full.png'))
    close(figFull)

    % zoomed on the middle third - reuses the same midpoint logic as the
    % zoomed velocity overlay in velocity_conditioning.m
    t_mid = t0 + duration/2;
    midWindow = [t_mid - midWindowSec/2, t_mid + midWindowSec/2];
    midWindow = [max(midWindow(1), t0), min(midWindow(2), t1)];

    figMid = newfig();
    plot(t_s, a, 'Color', [0.85 0.325 0.098], 'LineWidth', 0.8)
    xlim(midWindow)
    grid on
    box on
    xlabel('t (s)')
    ylabel('a (m/s^2)')
    title(sprintf('%s: Acceleration a(t), zoomed mid-trial (t = %.2f-%.2f s)', ...
        condLabel, midWindow(1), midWindow(2)))
    saveas(figMid, fullfile(outDir, 'acceleration_zoomed.png'))
    close(figMid)

    % start of trial (accel-from-rest transient)
    startWindow = [t0, min(t0 + edgeWindowSec, t1)];

    figStart = newfig();
    plot(t_s, a, 'Color', [0.85 0.325 0.098], 'LineWidth', 0.8)
    xlim(startWindow)
    grid on
    box on
    xlabel('t (s)')
    ylabel('a (m/s^2)')
    title(sprintf('%s: Acceleration a(t), start of trial (t = %.2f-%.2f s)', ...
        condLabel, startWindow(1), startWindow(2)))
    saveas(figStart, fullfile(outDir, 'acceleration_start.png'))
    close(figStart)

    % end of trial (decel-to-stop transient)
    endWindow = [max(t1 - edgeWindowSec, t0), t1];

    figEnd = newfig();
    plot(t_s, a, 'Color', [0.85 0.325 0.098], 'LineWidth', 0.8)
    xlim(endWindow)
    grid on
    box on
    xlabel('t (s)')
    ylabel('a (m/s^2)')
    title(sprintf('%s: Acceleration a(t), end of trial (t = %.2f-%.2f s)', ...
        condLabel, endWindow(1), endWindow(2)))
    saveas(figEnd, fullfile(outDir, 'acceleration_end.png'))
    close(figEnd)

    % raw vs filtered acceleration, same start/end windows - just a visual
    % check for filtfilt artifacts, judged by eye
    figStartCmp = newfig();
    plot(t_s, a_raw, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8)
    hold on
    plot(t_s, a, 'Color', [0.85 0.325 0.098], 'LineWidth', 1.2)
    xlim(startWindow)
    grid on
    box on
    xlabel('t (s)')
    ylabel('a (m/s^2)')
    title(sprintf('%s: a from raw v vs. a from filtered v, start of trial (t = %.2f-%.2f s)', ...
        condLabel, startWindow(1), startWindow(2)))
    legend({'a from raw v', 'a from filtered v'}, 'Box', 'off')
    saveas(figStartCmp, fullfile(outDir, 'acceleration_start_comparison.png'))
    close(figStartCmp)

    figEndCmp = newfig();
    plot(t_s, a_raw, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8)
    hold on
    plot(t_s, a, 'Color', [0.85 0.325 0.098], 'LineWidth', 1.2)
    xlim(endWindow)
    grid on
    box on
    xlabel('t (s)')
    ylabel('a (m/s^2)')
    title(sprintf('%s: a from raw v vs. a from filtered v, end of trial (t = %.2f-%.2f s)', ...
        condLabel, endWindow(1), endWindow(2)))
    legend({'a from raw v', 'a from filtered v'}, 'Box', 'off')
    saveas(figEndCmp, fullfile(outDir, 'acceleration_end_comparison.png'))
    close(figEndCmp)

    % peak magnitude + zero crossing count in each window, just a rough
    % feel for how oscillatory each signal looks, not a real verdict
    startMask = t_s >= startWindow(1) & t_s <= startWindow(2);
    endMask   = t_s >= endWindow(1)   & t_s <= endWindow(2);

    peakStart_raw  = max(abs(a_raw(startMask)));
    peakStart_filt = max(abs(a(startMask)));
    peakEnd_raw    = max(abs(a_raw(endMask)));
    peakEnd_filt   = max(abs(a(endMask)));

    zcStart_raw  = nnz(diff(sign(a_raw(startMask))) ~= 0);
    zcStart_filt = nnz(diff(sign(a(startMask))) ~= 0);
    zcEnd_raw    = nnz(diff(sign(a_raw(endMask))) ~= 0);
    zcEnd_filt   = nnz(diff(sign(a(endMask))) ~= 0);

    fprintf('%s: [start window] peak|a_raw| = %.4f m/s^2, peak|a_filt| = %.4f m/s^2, zero-crossings a_raw/a_filt = %d/%d\n', ...
        condName, peakStart_raw, peakStart_filt, zcStart_raw, zcStart_filt);
    fprintf('%s: [end window]   peak|a_raw| = %.4f m/s^2, peak|a_filt| = %.4f m/s^2, zero-crossings a_raw/a_filt = %d/%d\n', ...
        condName, peakEnd_raw, peakEnd_filt, zcEnd_raw, zcEnd_filt);

    fprintf('%s: done. Outputs in %s\n', condName, outDir);

    % rough steady-state check (a_x ~ 0), still provisional, see note above
    fprintf('%s: [rough steady-state window, not final phase segmentation]\n', condName);
    for w = 1:prov_nWin
        p = prov_trimFracs(w);
        prov_lo = t0 + p * duration;
        prov_hi = t1 - p * duration;
        prov_mask = (t_s >= prov_lo) & (t_s <= prov_hi);

        prov_mean_ax   = mean(a(prov_mask));
        prov_std_ax    = std(a(prov_mask));
        prov_maxabs_ax = max(abs(a(prov_mask)));

        prov_r = prov_r + 1;
        prov_trial{prov_r}     = condLabel;
        prov_windowPct(prov_r) = prov_midPct(w);
        prov_meanAx(prov_r)    = prov_mean_ax;
        prov_stdAx(prov_r)     = prov_std_ax;
        prov_maxAbsAx(prov_r)  = prov_maxabs_ax;

        fprintf(['%s:   middle %2.0f%% (t = %6.2f-%6.2f s): ', ...
            'mean a_x = %+.5f, std a_x = %.5f, max|a_x| = %.5f m/s^2\n'], ...
            condName, prov_midPct(w), prov_lo, prov_hi, ...
            prov_mean_ax, prov_std_ax, prov_maxabs_ax);
    end
end

% full-trial peak table across all conditions, both signals - the single
% largest magnitude here gets adopted as a_max, used symmetrically in
% both directions for the stability analysis
peakTable = table(peakCond, peakSignal, peakValue, ...
    'VariableNames', {'condition', 'signal', 'abs_a_max_mps2'});
writetable(peakTable, fullfile(resultsRoot, 'peak_accel_table.csv'))

[a_max_adopted, idxMax] = max(peakValue);
fprintf('\n=== Peak acceleration summary (all conditions) ===\n');
disp(peakTable)
fprintf('Adopted a_max = %.4f m/s^2 (%s, %s), used symmetrically in both directions.\n', ...
    a_max_adopted, peakCond{idxMax}, peakSignal{idxMax});

% rough steady-state window table, all conditions. columns:
%   trial            - condition label
%   window_pct       - middle fraction of trial duration kept (50/60/40%)
%   mean_ax_mps2     - mean a(t) in that window
%   std_ax_mps2      - std a(t) in that window
%   max_abs_ax_mps2  - max |a(t)| in that window (high value = window
%                      probably still has some transient in it)
%
% if mean_ax stays small across all three windows per trial, the a_x ~ 0
% assumption isn't too sensitive to where exactly the window edges fall.
provTable = table(prov_trial, prov_windowPct, prov_meanAx, prov_stdAx, prov_maxAbsAx, ...
    'VariableNames', {'trial', 'window_pct', 'mean_ax_mps2', 'std_ax_mps2', 'max_abs_ax_mps2'});
writetable(provTable, fullfile(resultsRoot, 'PROVISIONAL_steady_state_window_accel.csv'))

fprintf('\n=== rough steady-state window a_x (all conditions) ===\n');
disp(provTable)


function f = newfig()
% forces light theme so saved PNGs don't depend on the user's MATLAB desktop theme
    f = figure('Visible', 'off', 'Color', 'w');
    try
        f.Theme = 'light';   % R2025a+
    catch
        % older MATLAB, no theme support - white background above is enough
    end
end

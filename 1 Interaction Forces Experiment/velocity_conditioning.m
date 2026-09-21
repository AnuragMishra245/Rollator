clear
close all
clc

files = {
    'walker_straight_1_constvel-0.2ms_velocity.csv',  'const_v', 'Constant velocity 0.2 m/s';
    'walker_straight_2_adm-dv20-dw10_velocity.csv',   'adm1',    'Admittance 1 (dv=20, dw=10)';
    'walker_straight_3_adm-dv200-dw100_velocity.csv', 'adm2',    'Admittance 2 (dv=200, dw=100)'
    };

nCond = size(files, 1);
scriptDir = fileparts(mfilename('fullpath'));
dataRoot = fullfile(scriptDir, 'Interaction Force Measured Data');  % raw walker CSVs
resultsRoot = fullfile(scriptDir, 'Filtered Data');

cutoffCandidates = 1:20;   % Hz, Winter residual sweep
butterOrder = 2;
noise_region_hz = [10 20]; % Hz, high-cutoff range used for the Winter noise-floor fit
zoomWindowSec = 2;         % width of zoomed raw-vs-filtered overlay, s
transientWindowSec = 4;    % width of accel/decel overlay windows, s

summaryCond = cell(nCond, 1);
summaryCutoff = nan(nCond, 1);
summarySlopeA = nan(nCond, 1);
summaryInterceptB = nan(nCond, 1);

for k = 1:nCond
    srcFile  = files{k,1};
    condName = files{k,2};
    condLabel = files{k,3};

    outDir = fullfile(resultsRoot, condName);
    if ~exist(outDir, 'dir')
        mkdir(outDir)
    end

    fprintf('\n=== %s (%s) ===\n', condLabel, condName);

    % NaN/monotonic-time/range checks live in the data-QA notebook, not repeated here
    [t_s, vx, ~] = load_velocity_file(fullfile(dataRoot, srcFile));

    dt = diff(t_s);
    fs_mean = 1 / mean(dt);

    summaryFile = fullfile(outDir, 'sample_rate_summary.txt');
    fid = fopen(summaryFile, 'w');
    fprintf(fid, 'Condition: %s (%s)\n', condLabel, condName);
    fprintf(fid, 'Source file: %s\n', srcFile);
    fprintf(fid, 'N samples: %d\n', numel(t_s));
    fprintf(fid, 'Duration: %.3f s\n', t_s(end) - t_s(1));
    fprintf(fid, 'Mean fs: %.4f Hz\n', fs_mean);
    fclose(fid);

    % Winter cutoff selection
    residual = nan(size(cutoffCandidates));
    for ic = 1:numel(cutoffCandidates)
        fc = cutoffCandidates(ic);
        Wn = fc / (fs_mean/2);
        if Wn >= 1
            continue  % above Nyquist for this trial's fs
        end
        [b, a] = butter(butterOrder, Wn, 'low');
        vx_f = filtfilt(b, a, vx);
        residual(ic) = rms(vx - vx_f);
    end

    validIdx = ~isnan(residual);
    fc_valid = cutoffCandidates(validIdx);
    res_valid = residual(validIdx);

    % Winter (1990): fit a line to the high-cutoff part of the residual
    % curve (dominated by noise, roughly linear there), extrapolate to
    % cutoff = 0, and the intercept b is the noise-floor estimate. Chosen
    % cutoff is the lowest swept value where the actual residual drops to
    % or below b.
    noiseIdx = fc_valid >= noise_region_hz(1) & fc_valid <= noise_region_hz(2);
    if sum(noiseIdx) < 2
        error('velocity_conditioning:NoiseRegion', ...
            '%s: fewer than 2 valid sweep points in noise_region_hz [%g %g]', ...
            condName, noise_region_hz(1), noise_region_hz(2));
    end
    fitCoeffs = polyfit(fc_valid(noiseIdx), res_valid(noiseIdx), 1);
    a_slope = fitCoeffs(1);
    b_intercept = fitCoeffs(2);

    crossIdx = find(res_valid <= b_intercept, 1, 'first');
    if isempty(crossIdx)
        warning('%s: residual curve never dropped to the extrapolated noise floor; using highest swept cutoff', condName);
        crossIdx = numel(res_valid);
    end
    cutoff_chosen = fc_valid(crossIdx);

    fprintf('%s: chosen cutoff = %d Hz (Winter cutoff, regression method)\n', ...
        condName, cutoff_chosen);

    figResid = newfig();
    plot(fc_valid, res_valid, 'o-', 'LineWidth', 1.2)
    hold on
    fitLineX = [0, fc_valid(end)];
    fitLineY = a_slope*fitLineX + b_intercept;
    plot(fitLineX, fitLineY, 'k--', 'LineWidth', 1)
    plot(0, b_intercept, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 6)
    plot(cutoff_chosen, res_valid(crossIdx), 'r*', 'MarkerSize', 10)
    grid on
    box on
    xlabel('Cutoff frequency (Hz)')
    ylabel('Residual RMS (vx_{raw} - vx_{filt})')
    title(sprintf('%s: Residual vs. Cutoff (Winter method)', condLabel))
    legend({'residual', 'noise-region fit (extrapolated)', 'intercept b', ...
        sprintf('Winter cutoff = %d Hz', cutoff_chosen)}, 'Box', 'off', 'Location', 'best')
    saveas(figResid, fullfile(outDir, 'residual_curve.png'))
    close(figResid)

    residTable = table(fc_valid(:), res_valid(:), 'VariableNames', {'cutoff_Hz', 'residual_rms'});
    writetable(residTable, fullfile(outDir, 'residual_data.csv'))

    fid = fopen(summaryFile, 'a');
    fprintf(fid, '\nCutoff selection\n');
    fprintf(fid, 'Chosen cutoff (Winter cutoff): %d Hz\n', cutoff_chosen);
    fprintf(fid, 'Method: Winter (1990) residual analysis. Noise region = [%g %g] Hz.\n', ...
        noise_region_hz(1), noise_region_hz(2));
    fprintf(fid, 'Noise-region linear fit: residual = a*cutoff + b, a = %.6f, b = %.6f\n', a_slope, b_intercept);
    fprintf(fid, 'Selected cutoff = lowest swept cutoff where residual first drops to/below b (no interpolation).\n');
    fprintf(fid, 'Heuristic - check against residual_curve.png / residual_data.csv before locking in.\n');
    fclose(fid);

    summaryCond{k} = condName;
    summaryCutoff(k) = cutoff_chosen;
    summarySlopeA(k) = a_slope;
    summaryInterceptB(k) = b_intercept;

    % apply filter, save filtered velocity
    Wn_chosen = cutoff_chosen / (fs_mean/2);
    [b, a] = butter(butterOrder, Wn_chosen, 'low');
    vx_filt = filtfilt(b, a, vx);

    figOverlay = newfig();
    plot(t_s, vx, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8)
    hold on
    plot(t_s, vx_filt, 'Color', [0 0.447 0.741], 'LineWidth', 1.2)
    grid on
    box on
    xlabel('t (s)')
    ylabel('v_x (m/s)')
    title(sprintf('%s: Raw vs. Filtered v_x (f_c = %d Hz)', condLabel, cutoff_chosen))
    legend({'raw', 'filtered'}, 'Box', 'off')
    saveas(figOverlay, fullfile(outDir, 'filter_overlay.png'))
    close(figOverlay)

    % zoomed on the trial's temporal midpoint (no steady-state segmentation
    % exists elsewhere in the project yet, this just avoids the start/stop
    % transients) with a fixed width to inspect the filter over roughly a
    % gait cycle
    t_mid = t_s(1) + (t_s(end) - t_s(1)) / 2;
    zoomWindow = [t_mid - zoomWindowSec/2, t_mid + zoomWindowSec/2];

    figZoom = newfig();
    plot(t_s, vx, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8)
    hold on
    plot(t_s, vx_filt, 'Color', [0 0.447 0.741], 'LineWidth', 1.2)
    xlim(zoomWindow)
    grid on
    box on
    xlabel('t (s)')
    ylabel('v_x (m/s)')
    title(sprintf('%s: Raw vs. Filtered v_x, zoomed (t = %.1f-%.1f s, f_c = %d Hz)', ...
        condLabel, zoomWindow(1), zoomWindow(2), cutoff_chosen))
    legend({'raw', 'filtered'}, 'Box', 'off')
    saveas(figZoom, fullfile(outDir, 'filter_overlay_zoomed.png'))
    close(figZoom)

    % fixed windows at trial start/end. transientWindowSec (4 s) came from
    % eyeballing the raw vx ramps across the three conditions - roughly
    % 90% of steady-state by 1-3.2 s in, below 10% by 1-1.1 s from the end,
    % so 4 s gives some margin either way
    accelWindow = [t_s(1), t_s(1) + transientWindowSec];
    decelWindow = [t_s(end) - transientWindowSec, t_s(end)];

    figAccel = newfig();
    plot(t_s, vx, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8)
    hold on
    plot(t_s, vx_filt, 'Color', [0 0.447 0.741], 'LineWidth', 1.2)
    xlim(accelWindow)
    grid on
    box on
    xlabel('t (s)')
    ylabel('v_x (m/s)')
    title(sprintf('%s: Raw vs. Filtered v_x, acceleration phase (t = %.1f-%.1f s, f_c = %d Hz)', ...
        condLabel, accelWindow(1), accelWindow(2), cutoff_chosen))
    legend({'raw', 'filtered'}, 'Box', 'off')
    saveas(figAccel, fullfile(outDir, 'filter_overlay_accel.png'))
    close(figAccel)

    figDecel = newfig();
    plot(t_s, vx, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8)
    hold on
    plot(t_s, vx_filt, 'Color', [0 0.447 0.741], 'LineWidth', 1.2)
    xlim(decelWindow)
    grid on
    box on
    xlabel('t (s)')
    ylabel('v_x (m/s)')
    title(sprintf('%s: Raw vs. Filtered v_x, deceleration phase (t = %.1f-%.1f s, f_c = %d Hz)', ...
        condLabel, decelWindow(1), decelWindow(2), cutoff_chosen))
    legend({'raw', 'filtered'}, 'Box', 'off')
    saveas(figDecel, fullfile(outDir, 'filter_overlay_decel.png'))
    close(figDecel)

    outTable = table(t_s(:), vx(:), vx_filt(:), 'VariableNames', {'t_s', 'vx_raw', 'vx_filt'});
    writetable(outTable, fullfile(outDir, 'filtered_velocity.csv'))

    fprintf('%s: done. Outputs in %s\n', condName, outDir);
end

fprintf('\n=== Winter-cutoff summary (all conditions) ===\n');
fprintf('%-10s %12s %14s %14s\n', 'Condition', 'Cutoff (Hz)', 'Slope a', 'Intercept b');
for k = 1:nCond
    fprintf('%-10s %12d %14.6f %14.6f\n', ...
        summaryCond{k}, summaryCutoff(k), summarySlopeA(k), summaryInterceptB(k));
end


function f = newfig()
% forces light theme so saved PNGs don't depend on the user's MATLAB desktop theme
    f = figure('Visible', 'off', 'Color', 'w');
    try
        f.Theme = 'light';   % R2025a+
    catch
        % older MATLAB, no theme support - white background above is enough
    end
end

function [t_s, vx, wz] = load_velocity_file(filename)
% Reads a walker velocity CSV (t_s, vx, wz). Some exports come through as
% a normal 3-column CSV, others get packed into one text column (locale
% delimiter mismatch in Excel) - handles both.
    T = readtable(filename, 'Delimiter', ',');

    if width(T) >= 3
        t_s = T.(1);
        vx  = T.(2);
        wz  = T.(3);
    elseif width(T) == 1
        raw = string(T.(1));
        parts = split(raw, ',');
        t_s = str2double(parts(:,1));
        vx  = str2double(parts(:,2));
        wz  = str2double(parts(:,3));
    else
        error('load_velocity_file:UnexpectedFormat', ...
            '%s: unexpected column count (%d)', filename, width(T));
    end

    t_s = t_s(:);
    vx  = vx(:);
    wz  = wz(:);
end

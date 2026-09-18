%% Complete ECG Analysis: DF1 Algorithm & Statistical Validation
clc; clear; close all;

% =========================================================
% 1. SETUP & DATA IMPORT
% =========================================================
filename = 'part a.txt';
raw_data = readmatrix(filename); 

% Extraction: Assuming Lead II is Column 3
sig_lead2 = raw_data(:, 3); 
Fs = 500; 

% =========================================================
% 2. SIGNAL PRE-PROCESSING
% =========================================================
% Bandpass Filter (0.5 - 50 Hz) to remove noise/drift
bp_freqs = [0.5 50];
[num, den] = butter(2, bp_freqs / (Fs/2), 'bandpass');
sig_clean  = filtfilt(num, den, sig_lead2);

% --- DF1 Algorithm Steps ---
% Step A: Differentiation (8-sample delay)
h_diff = [1, zeros(1, 7), -1];
sig_diff = filter(h_diff, 1, sig_clean);

% Step B: Smoothing (Low Pass)
h_smooth = [1 4 6 4 1];
sig_integrated = filter(h_smooth, 1, sig_diff);

% =========================================================
% 3. R-PEAK DETECTION ENGINE
% =========================================================
% Thresholds & Windowing
TH_UP   = 1.7;
TH_DOWN = -1.7;
WIN_MS  = 160;
win_len = round((WIN_MS/1000) * Fs);

detected_peaks = [];
N = length(sig_integrated);
idx = 20; 

while idx <= N - win_len
    if sig_integrated(idx) > TH_UP
        
        win_end = min(N, idx + win_len - 1);
        search_window = idx : win_end;
        window_vals = sig_integrated(search_window);
        
        % Count Threshold Crossings
        cross_counter = 1; 
        state = 1; 
        
        for v = window_vals'
            if (state == 1 && v < TH_DOWN)
                cross_counter = cross_counter + 1;
                state = -1;
            elseif (state == -1 && v > TH_UP)
                cross_counter = cross_counter + 1;
                state = 1;
            end
        end
        
        % Validate QRS Complex
        if cross_counter >= 2 && cross_counter <= 4
            [~, local_max] = max(sig_clean(search_window));
            true_peak_loc = search_window(1) + local_max - 1;
            detected_peaks(end+1) = true_peak_loc; 
        end
        
        idx = win_end;
    else
        idx = idx + 1;
    end
end

% =========================================================
% 4. SEGMENTATION & STATISTICAL ANALYSIS
% =========================================================
t_axis = (0:length(sig_clean)-1) / Fs;

% Define Segments: [Start, End]
roi_defs = [0 30; 31 41; 42 52]; 
labels   = {'Seated','Standing', 'Deep Breathing'};

results = struct();

fprintf('\n--- Segment Analysis Results ---\n');

for k = 1:3
    t_start = roi_defs(k, 1);
    t_end   = roi_defs(k, 2);
    
    % Find peaks strictly within this time window
    peak_times = t_axis(detected_peaks);
    valid_indices = detected_peaks(peak_times >= t_start & peak_times <= t_end);
    
    % Calculate RR Intervals (in seconds)
    rr_intervals = diff(t_axis(valid_indices));
    
    % --- CORRECTION HERE: Artifact Rejection ---
    % Filter out impossible RR intervals (noise spikes)
    % We assume HR cannot be > 200 BPM (Interval < 0.3s)
    rr_clean = rr_intervals(rr_intervals > 0.3);
    
    % If segment is too noisy and all peaks removed, avoid crash
    if isempty(rr_clean)
        hr_bpm = 0;
    else
        hr_bpm = 60 ./ rr_clean;
    end
    % -------------------------------------------
    
    % Store Statistics
    results(k).name = labels{k};
    results(k).hr_mean = mean(hr_bpm);
    results(k).hr_std  = std(hr_bpm);
    results(k).bpm_data = hr_bpm; 
    results(k).num_cycles = length(hr_bpm);
    
    % Print Output
    fprintf('Segment %d (%s): HR = %.2f bpm, STD = %.2f bpm, Cycles = %d\n', ...
        k, labels{k}, results(k).hr_mean, results(k).hr_std, results(k).num_cycles);
end

% =========================================================
% 5. STATISTICAL TEST (Unpaired t-test)
% =========================================================
data_seated = results(2).bpm_data;
data_breathing = results(3).bpm_data;

[h, p_val, ci, stats] = ttest2(data_breathing, data_seated, 'Vartype', 'unequal');

fprintf('\n--- Statistical Comparison (Seated vs. Deep Breathing) ---\n');
fprintf('Method: Two-sample t-test (Welch''s approximation)\n');
fprintf('t-statistic = %.4f\n', stats.tstat);
fprintf('Degrees of Freedom (df) = %.2f\n', stats.df);
fprintf('P-value = %.8f\n', p_val);

if h == 1
    fprintf('Result: Significant difference detected (p < 0.05)\n');
else
    fprintf('Result: No significant difference found.\n');
end
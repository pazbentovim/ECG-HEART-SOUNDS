%% Advanced ECG Analysis: Segmentation, SNR & Ensemble Averaging
clc; clear; close all;

% --- 1. Data Ingestion & Preprocessing ---
raw_table = readmatrix('part b.txt');
sig_raw   = raw_table(:, 3); % Extract Lead II
Fs        = 500; 

% Bandpass Filter Setup (Butterworth, 2nd order, 0.5-50Hz)
freq_band = [0.5 50];
[num, den] = butter(2, freq_band / (Fs/2), 'bandpass');

% Apply Zero-Phase Filtering
sig_filt = filtfilt(num, den, sig_raw);

% Trim edges (1 second from start/end) to remove transient artifacts
trim_idx = Fs * 1;
main_signal = sig_filt(trim_idx+1 : end-trim_idx);
t_vec = (0:length(main_signal)-1) / Fs;

% --- 2. QRS Detection (DF1 Algorithm) ---
% Stage A: Derivative Filter
h_deriv = [1, zeros(1, 7), -1];
s_deriv = filter(h_deriv, 1, main_signal);

% Stage B: Smoothing Filter
h_smooth = [1 4 6 4 1];
s_integ  = filter(h_smooth, 1, s_deriv);

% Detection Logic
TH_UP   = 1.7; 
TH_DOWN = -1.7;
WIN_WIDTH_MS = 160;
win_len = round((WIN_WIDTH_MS/1000) * Fs);

detected_indices = [];
N = length(s_integ);
idx = 20; 

while idx <= N - win_len
    if s_integ(idx) > TH_UP
        % Define Analysis Window
        win_end = min(N, idx + win_len - 1);
        window_range = idx : win_end;
        
        % Check Crossing Pattern
        cross_count = 1; 
        state = 1; 
        
        vals = s_integ(window_range);
        for k = 2:length(vals)
            v = vals(k);
            if (state == 1 && v < TH_DOWN)
                cross_count = cross_count + 1;
                state = -1;
            elseif (state == -1 && v > TH_UP)
                cross_count = cross_count + 1;
                state = 1;
            end
        end
        
        % Valid QRS condition
        if cross_count >= 2 && cross_count <= 4
            [~, peak_offset] = max(main_signal(window_range));
            true_peak = window_range(1) + peak_offset - 1;
            detected_indices(end+1) = true_peak; 
        end
        
        idx = win_end;
    else
        idx = idx + 1;
    end
end

% --- 3. Statistical Filtering (Outlier Removal) ---
r_times_sec = t_vec(detected_indices);
rr_intervals = diff(r_times_sec);
bpm_inst = 60 ./ rr_intervals;

% Stats
mu_rr = mean(rr_intervals);
sigma_rr = std(rr_intervals);

% Identify outliers (> 2 standard deviations)
is_outlier = abs(rr_intervals - mu_rr) > (2 * sigma_rr);

% Clean Datasets
rr_clean = rr_intervals(~is_outlier);
bpm_clean = bpm_inst(~is_outlier);

fprintf('Total Detected Beats: %d\n', length(detected_indices));
fprintf('Rejected Outliers: %d\n', sum(is_outlier));
fprintf('Clean HR Mean: %.2f bpm\n', mean(bpm_clean));
fprintf('Clean HR Std:  %.2f bpm\n', std(bpm_clean));

% --- 4. Segmentation & Ensemble Averaging ---
% Window: -200ms to +300ms
pre_s = 0.2; post_s = 0.3;
pre_samp  = round(pre_s * Fs);
post_samp = round(post_s * Fs);

valid_beat_indices = find(~is_outlier) + 1;

beat_segments = [];
for k = 1:length(valid_beat_indices)
    beat_idx = valid_beat_indices(k);
    center_loc = detected_indices(beat_idx);
    
    idx_start = center_loc - pre_samp;
    idx_end   = center_loc + post_samp;
    
    if idx_start >= 1 && idx_end <= N
        beat_segments = [beat_segments; main_signal(idx_start:idx_end)']; 
    end
end

% Ensemble calculations
mean_template = mean(beat_segments, 1);
std_template  = std(beat_segments, 0, 1);
t_segment = (-pre_samp : post_samp) / Fs;

% Extract specific beats
idx_first = find(detected_indices - pre_samp >= 1, 1, 'first');
idx_last  = find(detected_indices + post_samp <= N, 1, 'last');

first_beat_sig = main_signal(detected_indices(idx_first)-pre_samp : detected_indices(idx_first)+post_samp)';
last_beat_sig  = main_signal(detected_indices(idx_last)-pre_samp  : detected_indices(idx_last)+post_samp)';

% --- 5. SNR Calculations ---
iso_dur_ms = 20;
iso_len = round((iso_dur_ms/1000) * Fs);
start_iso = pre_samp - round(0.05 * Fs);
range_iso = start_iso : (start_iso + iso_len - 1);

noise_avg   = std(mean_template(range_iso));
noise_first = std(first_beat_sig(range_iso));
% --- 5.5 Isoelectric Segment Statistics & R-Peak SNR ---

% Isoelectric segment definition (same logic as reference code)
iso_ms   = 20;                         
iso_samp = round((iso_ms/1000) * Fs);
iso_start = pre_samp - round(0.05 * Fs);
iso_end   = iso_start + iso_samp - 1;

% Safety check
if iso_start < 1
    error('Isoelectric segment starts before signal boundary');
end

% Extract isoelectric segments
iso_avg   = mean_template(iso_start:iso_end);
iso_first = first_beat_sig(iso_start:iso_end);
iso_last  = last_beat_sig(iso_start:iso_end);

% Noise statistics
std_avg   = std(iso_avg);
std_first = std(iso_first);
std_last  = std(iso_last);

% R-peak values
R_peak_avg   = max(mean_template);
R_peak_first = max(first_beat_sig);
R_peak_last  = max(last_beat_sig);

% SNR based on R-peak / noise STD
SNR_avg_R   = R_peak_avg   / std_avg;
SNR_first_R = R_peak_first / std_first;
SNR_last_R  = R_peak_last  / std_last;

% Print results
fprintf('\nIsoelectric segment statistics (%.1f ms):\n', iso_ms);
fprintf('Average beat: mean = %.4f mV, std = %.4f mV\n', mean(iso_avg), std_avg);
fprintf('First beat:   mean = %.4f mV, std = %.4f mV\n', mean(iso_first), std_first);
fprintf('Last beat:    mean = %.4f mV, std = %.4f mV\n', mean(iso_last), std_last);

fprintf('\nR-peak values and SNR (R-peak / noise STD):\n');
fprintf('Average beat: R-peak = %.4f mV, SNR = %.2f\n', R_peak_avg, SNR_avg_R);
fprintf('First beat:   R-peak = %.4f mV, SNR = %.2f\n', R_peak_first, SNR_first_R);
fprintf('Last beat:    R-peak = %.4f mV, SNR = %.2f\n', R_peak_last, SNR_last_R);

% Continuous SNR vectors
snr_vec_avg   = mean_template / noise_avg;
snr_vec_early = first_beat_sig / noise_first;
snr_gain      = snr_vec_avg - snr_vec_early;

% --- 6. Visualization with LEGENDS ---
% Color Palette
col_teal   = [0, 0.5, 0.5];
col_orange = [0.85, 0.325, 0.098];
col_purp   = [0.494, 0.184, 0.556];
col_grey   = [0.4, 0.4, 0.4];
col_blue   = [0, 0.447, 0.741];

% Figure A: SNR Comparisons
figure('Name', 'SNR Analysis', 'Color', 'w');
subplot(2,1,1);
plot(t_segment, snr_vec_avg, 'Color', col_teal, 'LineWidth', 1.5, 'DisplayName', 'Avg Beat SNR');
hold on;
plot(t_segment, snr_vec_early, 'Color', col_orange, 'LineWidth', 1.5, 'DisplayName', 'First Beat SNR');
title('Signal-to-Noise Ratio Comparison'); 
ylabel('SNR'); grid on; legend('show', 'Location', 'best');
xlim([min(t_segment) max(t_segment)]);

subplot(2,1,2);
plot(t_segment, snr_gain, 'Color', col_purp, 'LineWidth', 1.5, 'DisplayName', '\Delta SNR (Gain)');
title('SNR Improvement'); 
xlabel('Time [sec]'); ylabel('\Delta SNR'); 
grid on; legend('show', 'Location', 'best');
xlim([min(t_segment) max(t_segment)]);

% Figure B: Beat Morphology Comparison
figure('Name', 'Beat Morphology', 'Color', 'w');

% 1. Average Beat with STD bounds
subplot(3,1,1);
plot(t_segment, mean_template, 'Color', col_blue, 'LineWidth', 1.5, 'DisplayName', 'Mean Template'); 
hold on;
plot(t_segment, mean_template + std_template, ':', 'Color', col_orange, 'LineWidth', 1, 'DisplayName', '+1 STD');
plot(t_segment, mean_template - std_template, ':', 'Color', col_orange, 'LineWidth', 1, 'DisplayName', '-1 STD');
title('Ensemble Statistics');
ylabel('Amplitude [mV]'); grid on; axis tight;
legend('show', 'Location', 'northeast');

% 2. First vs Average
subplot(3,1,2);
plot(t_segment, first_beat_sig, 'Color', col_purp, 'LineWidth', 1.2, 'DisplayName', 'First Valid Beat'); 
hold on;
plot(t_segment, mean_template, '--', 'Color', col_grey, 'LineWidth', 1.2, 'DisplayName', 'Ensemble Mean');
title('Individual vs. Ensemble');
ylabel('Amplitude [mV]'); grid on; axis tight;
legend('show', 'Location', 'best');

% 3. Last vs Average
subplot(3,1,3);
plot(t_segment, last_beat_sig, 'Color', col_teal, 'LineWidth', 1.2, 'DisplayName', 'Last Valid Beat'); 
hold on;
plot(t_segment, mean_template, '--', 'Color', col_grey, 'LineWidth', 1.2, 'DisplayName', 'Ensemble Mean');
title('Individual vs. Ensemble');
xlabel('Time [sec]'); ylabel('Amplitude [mV]'); grid on; axis tight;
legend('show', 'Location', 'best');
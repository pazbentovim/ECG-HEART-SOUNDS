%% ECG R-Peak Detection (DF1 Implementation)
clc; clear; close all;

% --- 1. Data Ingestion ---
% Load the matrix and isolate Lead II (Column 3)
raw_input = readmatrix('part a.txt');
sig_raw   = raw_input(:, 3); 
Fs = 500; 

% --- 2. Pre-processing (Bandpass) ---
% 0.5-50Hz Butterworth Filter
params_bp = [0.5 50] / (Fs/2);
[num, den] = butter(2, params_bp, 'bandpass');

% Apply filter (Zero-phase)
sig_clean = filtfilt(num, den, sig_raw);

% --- 3. DF1 Transform Stages ---
% A. Derivative Stage (8-sample delay)
% Kernel: [1, 0, 0, 0, 0, 0, 0, 0, -1]
kernel_d  = [1, zeros(1, 7), -1];
s_deriv   = filter(kernel_d, 1, sig_clean);

% B. Integration/Smoothing Stage
% Kernel: Moving average approximation [1 4 6 4 1]
kernel_s  = [1 4 6 4 1];
s_integ   = filter(kernel_s, 1, s_deriv);

% --- 4. Peak Detection Engine ---
% Constants
UP_THRESH   = 1.7;
DOWN_THRESH = -1.7;
WIN_WIDTH   = 160; % ms
win_size    = round(WIN_WIDTH/1000 * Fs);

detected_indices = [];
N = length(s_integ);
curr_idx = 20; % Skip initial transient

while curr_idx <= N - win_size
    % Check for activation threshold
    if s_integ(curr_idx) > UP_THRESH
        
        % Define Analysis Window
        start_win = curr_idx;
        end_win   = min(N, start_win + win_size - 1);
        range_idxs = start_win:end_win;
        
        % Analyze Crossing Pattern (High -> Low -> High)
        % Using a for-loop instead of while for cleaner structure
        cross_count = 1; 
        polarity = 1; % 1 = Positive, -1 = Negative
        
        vals_in_window = s_integ(range_idxs);
        
        % Iterate through window to count flips
        for k = 2:length(vals_in_window)
            val = vals_in_window(k);
            if (polarity == 1 && val < DOWN_THRESH)
                cross_count = cross_count + 1;
                polarity = -1;
            elseif (polarity == -1 && val > UP_THRESH)
                cross_count = cross_count + 1;
                polarity = 1;
            end
        end
        
        % DF1 Validation Rule: 2 <= Crossings <= 4
        if cross_count >= 2 && cross_count <= 4
            % Find peak in the CLEAN signal (time alignment)
            [~, offset] = max(sig_clean(range_idxs));
            true_peak = range_idxs(1) + offset - 1;
            detected_indices(end+1) = true_peak; 
        end
        
        % Jump index to end of window
        curr_idx = end_win;
    else
        % Scan forward
        curr_idx = curr_idx + 1;
    end
end

% --- 5. Visualization (Matched to requirement) ---
t_axis = (0:length(sig_clean)-1) / Fs;

figure;
subplot(3,1,1); % Keeping the specific subplot structure requested

% Identical plotting commands to original specification
plot(t_axis, sig_clean); hold on;
plot(t_axis(detected_indices), sig_clean(detected_indices), ...
    'or', 'MarkerFaceColor', 'r', 'MarkerSize', 4);

xlabel('Time [sec]');
ylabel('Amplitude [mV]');
ylim([-0.5 1.2]);
title('Lead II – R peaks Detection');
legend('Filtered Lead II', 'Detected R peaks', 'Location', 'best');
grid on;
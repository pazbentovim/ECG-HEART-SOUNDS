%% DF1 Algorithm Implementation for R-Peak Detection
clc; clear; close all;

% --- 1. Signal Loading & Pre-processing ---
% Load data silently
raw_data = readmatrix('part a.txt'); 
% Extract Lead II (Column 3)
lead_ii_raw = raw_data(:, 3); 
Fs = 500; 

% Bandpass Filter (Butterworth, 2nd order, 0.5-50Hz)
% Used to remove baseline wander and muscle noise
bp_range = [0.5 50];
[num_bp, den_bp] = butter(2, bp_range/(Fs/2), 'bandpass');
sig_filtered = filtfilt(num_bp, den_bp, lead_ii_raw);

% --- 2. DF1 Transformations ---
% Step A: Derivative Filter (H(z) = 1 - z^-8)
% Enhance the slope of the QRS complex
h_deriv = [1, 0, 0, 0, 0, 0, 0, 0, -1]; 
sig_deriv = filter(h_deriv, 1, sig_filtered);

% Step B: Smoothing Filter (Moving Average)
% H(z) = (1 + z^-1 + z^-2 + z^-3 + z^-4)^2 (approx) -> [1 4 6 4 1]
h_integ = [1 4 6 4 1];
sig_integrated = filter(h_integ, 1, sig_deriv);

% --- 3. Global R-Peak Detection Logic ---
% Constants based on DF1
THRESH_UP = 1.7;
THRESH_DOWN = -1.7;
WIN_SIZE_MS = 160;
win_samples = round((WIN_SIZE_MS/1000) * Fs);

detected_indices = [];
N = length(sig_integrated);
idx = 10; % Safety margin for filter transient

% Scan the entire processed signal
while idx < (N - win_samples)
    if sig_integrated(idx) > THRESH_UP
        % Candidate block found
        search_window = idx : (idx + win_samples - 1);
        
        % Analyze crossings logic (Pattern: High -> Low -> High)
        % We check how many times the signal oscillates between thresholds
        segment_vals = sig_integrated(search_window);
        
        % Count transitions
        cross_count = 1; % Start with 1 (the initial breach)
        state = 1; % 1 = above positive, -1 = below negative
        
        for v = segment_vals'
            if (state == 1 && v < THRESH_DOWN)
                cross_count = cross_count + 1;
                state = -1;
            elseif (state == -1 && v > THRESH_UP)
                cross_count = cross_count + 1;
                state = 1;
            end
        end
        
        % Validate QRS (Must have 2-4 crossings)
        if cross_count >= 2 && cross_count <= 4
            % Find the peak location in the filtered signal (time alignment)
            [~, local_peak_offset] = max(sig_filtered(search_window));
            peak_loc = search_window(1) + local_peak_offset - 1;
            detected_indices(end+1) = peak_loc; 
        end
        
        % Jump window
        idx = idx + win_samples;
    else
        idx = idx + 1;
    end
end

% --- 4. Visualization with Specific Titles ---
% Time vector
t_full = (0:N-1) / Fs;

% Defined Segments and Custom Titles
% Format: [Start_Time, End_Time]
roi_defs = [0, 30;    % Segment 1
            31, 41;   % Segment 2
            42, 52];  % Segment 3

% EXACT requested titles
custom_titles = {
    'Lead II - R peaks Detection (Seated)', ...
    'Lead II - R peaks Detection (Standing)', ...
    'Lead II - R peaks Detection (Deep Breathing)'
};

figure('Name', 'DF1 Detection Results', 'Color', 'w');

for k = 1:3
    subplot(3, 1, k);
    
    % Get window limits
    t1 = roi_defs(k, 1);
    t2 = roi_defs(k, 2);
    
    % Plot Signal (Filtered)
    plot(t_full, sig_filtered, 'LineWidth', 0.8); 
    hold on;
    
    % Overlay Detected Peaks (Only those within current view)
    valid_peaks = detected_indices(t_full(detected_indices) >= t1 & t_full(detected_indices) <= t2);
    plot(t_full(valid_peaks), sig_filtered(valid_peaks), ...
         'or', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    
    % Styling
    xlim([t1, t2]);
    ylim([-0.5 1.2]); % Uniform Y-axis
    grid on;
    
    title(custom_titles{k});
    ylabel('Amplitude [mV]');
    
    % Legend only on first plot
    if k == 1
        legend('ECG (Filtered)', 'Detected R-peaks', 'Location', 'best');
    end
    
    % X-label only on bottom plot
    if k == 3
        xlabel('Time [sec]');
    end
end
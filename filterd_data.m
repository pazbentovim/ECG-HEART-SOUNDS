%% Comprehensive ECG Report: Full Analysis & ROI Zoom
clc; clear; close all;

% =========================================================
% 1. DATA IMPORT & CONFIGURATION
% =========================================================
filename = 'part a.txt';
ecg_raw = readmatrix(filename); % Semicolon suppresses output
Fs = 500; 

% Channel Mapping: [Lead I (Col1), Lead III (Col2), Lead II (Col3)]
% Desired Plot Order: Lead I, Lead II, Lead III
lead_map = [1, 3, 2]; 
lead_labels = {'Lead I', 'Lead II', 'Lead III'};

% =========================================================
% 2. SIGNAL PROCESSING (FILTERING)
% =========================================================
% Design: 2nd Order Butterworth Bandpass (0.6 - 35 Hz)
f_pass = [0.6 35];
[b, a] = butter(2, f_pass / (Fs/2), 'bandpass');

% Apply Zero-Phase Filter
ecg_clean = filtfilt(b, a, ecg_raw);

% Global Time Vector
N = size(ecg_raw, 1);
t_vec = (0:N-1) / Fs;

% =========================================================
% 3. VISUALIZATION PART A: FULL SIGNAL (TIME DOMAIN)
% =========================================================
figure('Name', '1. Full Time Domain - Filtered', 'Color', 'w');
for k = 1:3
    subplot(3, 1, k);
    idx = lead_map(k); 
    
    plot(t_vec, ecg_clean(:, idx), 'LineWidth', 1);
    grid on;
    title(['Filtered Signal – ' lead_labels{k}]);
    ylabel('Amp [mV]');
    if k == 3, xlabel('Time [sec]'); end
    axis tight; 
    
    % Optional limit for Lead III stability
    if k == 3
        ylim([-0.2 0.6]);
    end
end
linkaxes(findall(gcf,'type','axes'), 'x');

% =========================================================
% 4. VISUALIZATION PART B: FREQUENCY DOMAIN (FFT)
% =========================================================
f_axis = (-ceil((N-1)/2) : floor((N-1)/2)) * (Fs/N); 
zoom_freq = 50; 
idx_roi_freq = f_axis >= -zoom_freq & f_axis <= zoom_freq;

data_sources = {ecg_raw, ecg_clean};
titles_fft = {'Frequency Domain - Raw Data', ...
              'Frequency Domain - Filtered Data'};

for d = 1:2
    current_data = data_sources{d};
    figure('Name', titles_fft{d}, 'Color', 'w');
    
    % קובע תיאור מצב הסיגנל (לפני / אחרי סינון)
    if d == 1
        signal_state = 'Before Filtering';
    else
        signal_state = 'After Filtering';
    end
    
    for k = 1:3
        subplot(3, 1, k);
        col_idx = lead_map(k);
        
        % Remove DC -> FFT -> Shift -> Normalize
        sig_ac = current_data(:, col_idx) - mean(current_data(:, col_idx)); 
        spec = fftshift(fft(sig_ac));
        mag_norm = abs(spec) / max(abs(spec));
        
        plot(f_axis(idx_roi_freq), mag_norm(idx_roi_freq), 'LineWidth', 1.2);
        grid on;
        title([lead_labels{k} ' – FFT (' signal_state ')']);
        ylabel('Norm. Mag');
        if k == 3, xlabel('Frequency [Hz]'); end
    end
end

% =========================================================
% 5. VISUALIZATION PART C: ROI ZOOM
% =========================================================
t_start = 10.2;
t_end   = 12.2;

idx_start = round(t_start * Fs) + 1;
idx_end   = round(t_end * Fs);
indices_roi = idx_start:idx_end;

t_segment = t_vec(indices_roi);

figure('Name', '4. ROI Zoom (3 Beats)', 'Color', 'w');
for k = 1:3
    subplot(3, 1, k);
    
    full_col = ecg_clean(:, lead_map(k));
    segment_data = full_col(indices_roi);
    
    plot(t_segment, segment_data, 'LineWidth', 1.2);
    grid on;
    title(sprintf('Filtered ECG – %s (3 beats)', lead_labels{k}));
    ylabel('Amplitude [mV]');
    
    if k == 3
        xlabel('Time [sec]');
    end
    
    xlim([t_start, t_end]);
end

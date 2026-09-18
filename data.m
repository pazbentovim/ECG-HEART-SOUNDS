%% Comprehensive ECG Report: Full Signal, Spectrum & ROI Analysis
clc; clear; close all;

% --- 1. Configuration & Data Import ---
filename = 'part a.txt';
raw_data = readmatrix(filename); 
Fs = 500; 

% Lead Mapping: [Lead I (Col1), Lead III (Col2), Lead II (Col3)]
% Desired Order: I, II, III -> Indices: [1, 3, 2]
channel_map = [1, 3, 2]; 
lead_titles = {'Lead I', 'Lead II', 'Lead III'};

% --- 2. Global Time & Frequency Calculations ---
N = size(raw_data, 1);
t_full = (0:N-1) / Fs;

% FFT Axis Setup (Robust for any N)
f_axis = (-ceil((N-1)/2) : floor((N-1)/2)) * (Fs/N);
zoom_freq = 60; % Focus FFT view on relevant range (+/- 60Hz)
idx_freq_roi = f_axis >= -zoom_freq & f_axis <= zoom_freq;

% --- 3. Region of Interest (ROI) Extraction ---
% Window: 10.2s to 12.2s
roi_times = [10.2, 12.2];
idx_start = round(roi_times(1) * Fs) + 1;
idx_end   = round(roi_times(2) * Fs);
roi_indices = idx_start:idx_end;

t_roi = t_full(roi_indices); % Time vector for the specific window

% --- 4. Visualization Engine ---
% Initialize 3 separate figures for clean reporting
hFig_Full = figure('Name', 'Full Time Domain', 'Color', 'w');
hFig_Freq = figure('Name', 'Frequency Domain (FFT)', 'Color', 'w');
hFig_Zoom = figure('Name', 'ROI Zoom (3 Beats)', 'Color', 'w');

for k = 1:3
    % Extract the specific lead data
    current_channel = raw_data(:, channel_map(k));
    
    % --- Plot A: Full Time Domain ---
    figure(hFig_Full);
    subplot(3, 1, k);
    plot(t_full, current_channel, 'LineWidth', 0.8);
    title([lead_titles{k} '  – Raw Data ECG']);
    ylabel('Amp [mV]');
    if k==3, xlabel('Time [sec]'); end
    axis tight; grid on;

    % --- Plot B: Frequency Domain (FFT) ---
    figure(hFig_Freq);
    subplot(3, 1, k);
    
    % Process FFT (Remove DC -> FFT -> Shift -> Normalize)
    sig_no_dc = current_channel - mean(current_channel);
    spec = fftshift(fft(sig_no_dc));
    mag_norm = abs(spec) / max(abs(spec));
    
    plot(f_axis(idx_freq_roi), mag_norm(idx_freq_roi) );
    title([lead_titles{k} ' – FFT of Raw Data']);
    ylabel('Norm. Mag');
    if k==3, xlabel('Frequency [Hz]'); end
    grid on;

    figure(hFig_Zoom);
    subplot(3, 1, k);
    
    % Extract just the segment
    segment_data = current_channel(roi_indices);
    
    plot(t_roi, segment_data,  'LineWidth', 1.2); 
    title([lead_titles{k} ' – Raw Data ECG (3 beats)']);
    ylabel('Amp [mV]');
    if k==3, xlabel('Time [sec]'); end
    grid on; xlim(roi_times);
end

% --- 5. Final Adjustments ---
% Sync axes for better interaction
figure(hFig_Full); linkaxes(findall(gcf, 'type', 'axes'), 'x');
figure(hFig_Zoom); linkaxes(findall(gcf, 'type', 'axes'), 'x');
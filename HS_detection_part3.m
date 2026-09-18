%% Electro-Mechanical Synchronization Analysis (ECG & PCG)
% Updated Visuals: Purple/Teal Theme
clc; clear; close all;

% =========================================================
% 1. Initialization & Data Loading
% =========================================================
raw_data = readmatrix('part c.txt');
sig_pcg = raw_data(:, 1); 
sig_ecg = raw_data(:, 2); 
Fs = 500; 

% Create Time Vector
N = length(sig_pcg);
time_vec = (0:N-1) / Fs;

% Define Segments
Segments(1).Name = 'Resting State';
Segments(1).Idx  = time_vec >= 0 & time_vec <= 20;

Segments(2).Name = 'Post-Exercise';
Segments(2).Idx  = time_vec >= 22 & time_vec <= 32;

% =========================================================
% 2. Signal Pre-processing
% =========================================================
[b_ecg, a_ecg] = butter(2, [0.5 50]/(Fs/2), 'bandpass');
[b_pcg, a_pcg] = butter(2, [25 200]/(Fs/2), 'bandpass');

% --- Main Processing Loop ---
results = [];

for k = 1:2
    % Extract Segment
    mask = Segments(k).Idx;
    t_seg = time_vec(mask);
    
    % Apply Filters
    ecg_filt = filtfilt(b_ecg, a_ecg, sig_ecg(mask));
    pcg_filt = filtfilt(b_pcg, a_pcg, sig_pcg(mask));
    
    % Normalize PCG (-1 to 1)
    pcg_norm = pcg_filt / max(abs(pcg_filt));
    
    % --- Feature Extraction ---
    % 1. Detect R-Peaks
    [r_locs, r_amps] = find_r_peaks(ecg_filt, Fs);
    
    % 2. Detect Heart Sounds (Gating Method)
    [s1_locs, s2_locs] = find_heart_sounds(pcg_norm, r_locs, Fs);
    
    % --- Metrics Calculation ---
    rr_intervals = diff(r_locs) / Fs;
    bpm = 60 ./ rr_intervals;
    
    num_cycles = min([length(r_locs), length(s1_locs), length(s2_locs)]) - 1;
    dt_r_s1 = zeros(1, num_cycles);
    dt_r_s2 = zeros(1, num_cycles);
    sys_dur = zeros(1, num_cycles); 
    dia_dur = zeros(1, num_cycles); 
    
    for i = 1:num_cycles
        dt_r_s1(i) = (s1_locs(i) - r_locs(i)) / Fs;
        dt_r_s2(i) = (s2_locs(i) - r_locs(i)) / Fs;
        sys_dur(i) = (s2_locs(i) - s1_locs(i)) / Fs;
        if i < num_cycles
            dia_dur(i) = (s1_locs(i+1) - s2_locs(i)) / Fs;
        end
    end
    dia_dur = dia_dur(1:end-1); 
    
    % Store Data
    Segments(k).Data.t = t_seg;
    Segments(k).Data.ecg = ecg_filt;
    Segments(k).Data.pcg = pcg_norm;
    Segments(k).Data.r = r_locs;
    Segments(k).Data.s1 = s1_locs;
    Segments(k).Data.s2 = s2_locs;
    
    Segments(k).Stats.Mean = [mean(bpm), mean(dt_r_s1), mean(dt_r_s2), mean(sys_dur), mean(dia_dur)];
    Segments(k).Stats.Std  = [std(bpm), std(dt_r_s1), std(dt_r_s2), std(sys_dur), std(dia_dur)];
end

% =========================================================
% 3. Results Table Generation
% =========================================================
vals_rest = Segments(1).Stats.Mean;
std_rest  = Segments(1).Stats.Std;
vals_act  = Segments(2).Stats.Mean;
std_act   = Segments(2).Stats.Std;
pct_change = ((vals_act - vals_rest) ./ vals_rest) * 100;

Metrics = {'Heart Rate [BPM]'; 'R-to-S1 Delay [sec]'; 'R-to-S2 Delay [sec]'; 'Systole (S1-S2) [sec]'; 'Diastole (S2-S1) [sec]'};
Str_Rest = compose('%.4f ± %.4f', vals_rest', std_rest');
Str_Act  = compose('%.4f ± %.4f', vals_act', std_act');
Str_Chg  = compose('%.2f %%', pct_change');

T = table(Str_Rest, Str_Act, Str_Chg, 'RowNames', Metrics, 'VariableNames', {'At_Rest', 'Post_Exercise', 'Change'});
disp(' '); disp('--- Electro-Mechanical Analysis Results ---'); disp(T);

% =========================================================
% 4. Visualization (New Color Scheme)
% =========================================================
% New Palette:
c_ecg = [0.494 0.184 0.556]; % Deep Purple
c_pcg = [0 0.5 0.5];         % Teal / Sea Green

for k = 1:2
    D = Segments(k).Data;
    figure('Color', 'w', 'Name', Segments(k).Name);
    t = linspace(D.t(1), D.t(end), length(D.ecg));
    
    % Subplot 1: Heart Sounds (PCG)
    subplot(2,1,1);
    plot(t, D.pcg, 'Color', c_pcg, 'LineWidth', 1.2); hold on;
    
    % Markers: Black for S1, Red for S2 (High Contrast)
    scatter(t(D.s1), D.pcg(D.s1)+0.08, 60, 'k', 'filled', '^');
    scatter(t(D.s2), D.pcg(D.s2)+0.08, 60, 'r', 'filled', 'v');
    
    title([Segments(k).Name ' - Phonocardiogram']);
    ylabel('Amplitude'); legend('PCG Signal', 'S1 Sound', 'S2 Sound', 'Location','northeast');
    grid on; grid minor; axis tight; ylim([-1.1 1.3]);
    
    % Subplot 2: ECG
    subplot(2,1,2);
    plot(t, D.ecg, 'Color', c_ecg, 'LineWidth', 1.2); hold on;
    scatter(t(D.r), D.ecg(D.r), 45, 'r', 'filled', 'o');
    
    % Visual Aid lines
    y_limits = ylim;
    for m = 1:min(6, length(D.s1)) % Show first 6 beats lines only
        xline(t(D.s1(m)), ':', 'Color', [0 0 0 0.4], 'LineWidth', 1);
    end
    
    title([Segments(k).Name ' - ECG Lead II']);
    xlabel('Time [sec]'); ylabel('mV');
    grid on; grid minor; axis tight;
end

%% --- Helper Functions ---

function [r_locs, r_amps] = find_r_peaks(ecg, fs)
    d_ecg = diff([0; ecg]);
    sq_ecg = d_ecg .^ 2;
    win_size = round(0.150 * fs); 
    int_ecg = movmean(sq_ecg, win_size);
    [~, locs_temp] = findpeaks(int_ecg, 'MinPeakHeight', mean(int_ecg)*2, 'MinPeakDistance', 0.3*fs);
    
    r_locs = zeros(size(locs_temp));
    r_amps = zeros(size(locs_temp));
    search_w = round(0.05*fs);
    for i = 1:length(locs_temp)
        idx = locs_temp(i);
        range = max(1, idx-search_w) : min(length(ecg), idx+search_w);
        [r_amps(i), rel_idx] = max(ecg(range));
        r_locs(i) = range(1) + rel_idx - 1;
    end
end

function [s1_idxs, s2_idxs] = find_heart_sounds(pcg, r_locs, fs)
    % Physiological Gating
    se = -1 * (pcg.^2) .* log(pcg.^2 + eps);
    env = movmean(se, 0.02 * fs); 
    env = env / max(env);
    
    s1_idxs = []; s2_idxs = [];
    
    for i = 1:length(r_locs)-1
        r_t = r_locs(i);
        next_r = r_locs(i+1);
        rr_dist = next_r - r_t;
        
        % S1 Window (20-120ms post R)
        w1_start = r_t + round(0.02 * fs);
        w1_end   = r_t + round(0.12 * fs);
        if w1_end > length(env), continue; end
        
        [~, mx1] = max(env(w1_start:w1_end));
        s1_cand = w1_start + mx1 - 1;
        
        if env(s1_cand) > 0.05
            s1_idxs(end+1) = s1_cand;
            
            % S2 Window (After S1 to 80% RR)
            w2_start = s1_cand + round(0.15 * fs); 
            w2_end   = r_t + round(0.75 * rr_dist);
            
            if w2_end <= length(env) && w2_start < w2_end
                [~, mx2] = max(env(w2_start:w2_end));
                s2_cand = w2_start + mx2 - 1;
                if env(s2_cand) > 0.05
                    s2_idxs(end+1) = s2_cand;
                end
            end
        end
    end
end
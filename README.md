[readme_md.md](https://github.com/user-attachments/files/32377467/readme_md.md)
# Quantitative Electrocardiogram (ECG) & Phonocardiogram (PCG) Analysis

This repository contains the digital signal processing (DSP) codebase, raw physiological data, and comprehensive research documentation for an advanced Biomedical Engineering study conducted at Tel Aviv University. The project investigates the critical synchronization between the electrical excitation and mechanical contraction of the human heart through the simultaneous acquisition and processing of ECG and PCG signals. 

Designed to extract reliable clinical metrics from noisy physiological recordings, this project demonstrates proficiency in biosignal processing, algorithm implementation, and statistical validation.

## Clinical & Engineering Highlights

* **Advanced Signal Pre-processing:** Implemented custom linear-phase FIR and 2nd-order zero-phase Butterworth bandpass filters to eliminate baseline wander (respiration artifacts) and attenuate high-frequency electromyographic (EMG) noise.
* **Automated QRS Detection Engine:** Developed a robust MATLAB implementation of the DF1 algorithm, utilizing signal differentiation, low-pass smoothing, and dynamic thresholding to identify R-peaks across varying physiological states.
* **Electro-Mechanical Synchronization:** Engineered a physiological gating algorithm using third-order Shannon Energy to accurately isolate and identify the first (S1) and second (S2) heart sounds from noisy acoustic data.
* **Statistical Rigor:** Applied ensemble averaging techniques to quantify Signal-to-Noise Ratio (SNR) improvements and utilized Welch's two-sample t-tests to validate the statistical significance of heart rate variability across different breathing paradigms.

## Repository Architecture

* `ECG & Heart Sounds.pdf`
  The complete, highly detailed final report containing physiological background, mathematical methodologies, graphical representations of the time/frequency domains, and extended clinical discussions.
* `filterd_data.m`
  MATLAB script responsible for primary signal cleaning. It applies a 0.6–35 Hz bandpass filter to raw ECG data and generates frequency domain spectrums via Fast Fourier Transform (FFT).
* `HR_T_TEST.m`
  Script executing the DF1 algorithm for feature extraction. It computes instantaneous beat-to-beat (RR) intervals, applies physiological artifact rejection (filtering HR > 200 BPM), and calculates statistical variances between resting and deep breathing states.
* `HS_detection_part3.m`
  The core script for electro-mechanical analysis. It normalizes PCG signals, applies the Shannon Energy envelope for S1/S2 detection, and maps mechanical acoustic events to electrical R-peaks to calculate systolic/diastolic durations.
* `part a.txt`, `part b.txt`, `part c.txt`
  Raw data matrices containing multi-lead ECG (Leads I, II, III) and stethoscope audio streams acquired at a 500 Hz sampling rate via a BIOPAC system.

## Key Experimental Results

The application of ensemble averaging successfully increased the ECG Signal-to-Noise Ratio (SNR) from 55.43 in a single raw beat to 67.71 in the cleaned ensemble. Furthermore, the synchronization algorithms successfully mapped the cardiac response to physical exertion, proving a strong electro-mechanical coupling under stress.

| Physiological Metric | At Rest (Mean ± STD) | Post-Exercise (Mean ± STD) | Shift |
| :--- | :--- | :--- | :--- |
| **Heart Rate** | 65.16 ± 3.20 BPM | 97.11 ± 5.69 BPM | +49.05% |
| **R-to-S1 Delay** | 0.0621 ± 0.0127 sec | 0.0341 ± 0.0068 sec | -45.03% |
| **Systole (S1 to S2)** | 0.3512 ± 0.0626 sec | 0.2587 ± 0.0497 sec | -26.34% |

## How to Run

1. Clone the repository and ensure MATLAB with the Signal Processing Toolbox is installed.
2. Ensure the raw data files (`part a.txt`, `part b.txt`, `part c.txt`) are located in the same working directory as the `.m` scripts.
3. Run `filterd_data.m` to observe the baseline correction and FFT transformations.
4. Run `HS_detection_part3.m` to generate the overlaid ECG/PCG synchronization plots and output the electro-mechanical timing matrices.

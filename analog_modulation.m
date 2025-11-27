% Analog Communication: AM Modulation Project
% DSB-SC and SSB-SC Modulation and Demodulation

clear all;
close all;
clc;

%% Parameters
fm = 1000;              % Message signal frequency (Hz)
fc = 10000;             % Carrier frequency (Hz)
fs = 100000;            % Sampling frequency (Hz)
t = 0:1/fs:0.01;        % Time vector (10ms)

%% Message Signal
m_t = cos(2*pi*fm*t);   % Message signal

%% Carrier Signal
c_t = cos(2*pi*fc*t);   % Carrier signal

%% DSB-SC Modulation
dsb_sc = m_t .* c_t;

%% SSB-SC Modulation (Upper Sideband)
% Using Hilbert transform method
m_t_hilbert = imag(hilbert(m_t));  % Hilbert transform of message
c_t_hilbert = sin(2*pi*fc*t);      % 90-degree phase shifted carrier

% SSB-SC (USB) = m(t)*cos(wc*t) - m_hat(t)*sin(wc*t)
ssb_sc_usb = m_t .* c_t - m_t_hilbert .* c_t_hilbert;

% SSB-SC (LSB) = m(t)*cos(wc*t) + m_hat(t)*sin(wc*t)
ssb_sc_lsb = m_t .* c_t + m_t_hilbert .* c_t_hilbert;

%% Demodulation - Coherent Detection

% DSB-SC Demodulation
dsb_demod = dsb_sc .* c_t;          % Multiply by carrier
[b, a] = butter(5, 2*fm/fs);        % Low-pass filter design
dsb_demod_filtered = 2*filter(b, a, dsb_demod);

% SSB-SC USB Demodulation
ssb_usb_demod = ssb_sc_usb .* c_t;
ssb_usb_demod_filtered = 2*filter(b, a, ssb_usb_demod);

% SSB-SC LSB Demodulation
ssb_lsb_demod = ssb_sc_lsb .* c_t;
ssb_lsb_demod_filtered = 2*filter(b, a, ssb_lsb_demod);

%% Frequency Domain Analysis
N = length(t);
f = (-N/2:N/2-1)*(fs/N);

% FFT of signals
M_f = fftshift(fft(m_t))/N;
DSB_f = fftshift(fft(dsb_sc))/N;
SSB_USB_f = fftshift(fft(ssb_sc_usb))/N;
SSB_LSB_f = fftshift(fft(ssb_sc_lsb))/N;

%% Plotting
figure('Position', [100, 100, 1200, 800]);

% Time Domain Signals
subplot(4,3,1);
plot(t*1000, m_t, 'b', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('Amplitude');
title('Message Signal m(t)');
grid on;

subplot(4,3,2);
plot(t*1000, c_t, 'r', 'LineWidth', 1);
xlabel('Time (ms)'); ylabel('Amplitude');
title('Carrier Signal c(t)');
grid on;
xlim([0 1]);

% DSB-SC
subplot(4,3,4);
plot(t*1000, dsb_sc, 'g', 'LineWidth', 1);
xlabel('Time (ms)'); ylabel('Amplitude');
title('DSB-SC Modulated Signal');
grid on;

subplot(4,3,5);
plot(f/1000, abs(M_f), 'b', 'LineWidth', 1.5);
xlabel('Frequency (kHz)'); ylabel('Magnitude');
title('Message Signal Spectrum');
grid on;
xlim([-fc/1000*1.5 fc/1000*1.5]);

subplot(4,3,6);
plot(f/1000, abs(DSB_f), 'g', 'LineWidth', 1.5);
xlabel('Frequency (kHz)'); ylabel('Magnitude');
title('DSB-SC Spectrum');
grid on;
xlim([-fc/1000*1.5 fc/1000*1.5]);

% SSB-SC USB
subplot(4,3,7);
plot(t*1000, ssb_sc_usb, 'm', 'LineWidth', 1);
xlabel('Time (ms)'); ylabel('Amplitude');
title('SSB-SC (USB) Modulated Signal');
grid on;

subplot(4,3,8);
plot(f/1000, abs(SSB_USB_f), 'm', 'LineWidth', 1.5);
xlabel('Frequency (kHz)'); ylabel('Magnitude');
title('SSB-SC (USB) Spectrum');
grid on;
xlim([-fc/1000*1.5 fc/1000*1.5]);

% SSB-SC LSB
subplot(4,3,9);
plot(t*1000, ssb_sc_lsb, 'c', 'LineWidth', 1);
xlabel('Time (ms)'); ylabel('Amplitude');
title('SSB-SC (LSB) Modulated Signal');
grid on;

subplot(4,3,10);
plot(f/1000, abs(SSB_LSB_f), 'c', 'LineWidth', 1.5);
xlabel('Frequency (kHz)'); ylabel('Magnitude');
title('SSB-SC (LSB) Spectrum');
grid on;
xlim([-fc/1000*1.5 fc/1000*1.5]);

% Demodulated Signals
subplot(4,3,11);
plot(t*1000, m_t, 'b', 'LineWidth', 1.5); hold on;
plot(t*1000, dsb_demod_filtered, 'g--', 'LineWidth', 1.5);
plot(t*1000, ssb_usb_demod_filtered, 'm:', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('Amplitude');
title('Demodulated Signals Comparison');
legend('Original', 'DSB-SC', 'SSB-USB');
grid on;

subplot(4,3,12);
plot(t*1000, m_t, 'b', 'LineWidth', 1.5); hold on;
plot(t*1000, ssb_lsb_demod_filtered, 'c--', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('Amplitude');
title('SSB-LSB Demodulated Signal');
legend('Original', 'SSB-LSB');
grid on;

sgtitle('Analog Communication: AM Modulation (DSB-SC and SSB-SC)', 'FontSize', 14, 'FontWeight', 'bold');

%% Performance Metrics
fprintf('=== Performance Metrics ===\n');
fprintf('DSB-SC Signal Power: %.4f\n', mean(dsb_sc.^2));
fprintf('SSB-SC USB Signal Power: %.4f\n', mean(ssb_sc_usb.^2));
fprintf('SSB-SC LSB Signal Power: %.4f\n', mean(ssb_sc_lsb.^2));
fprintf('\n');

% Calculate correlation with original signal
corr_dsb = corrcoef(m_t(100:end), dsb_demod_filtered(100:end));
corr_ssb_usb = corrcoef(m_t(100:end), ssb_usb_demod_filtered(100:end));
corr_ssb_lsb = corrcoef(m_t(100:end), ssb_lsb_demod_filtered(100:end));

fprintf('Demodulation Correlation:\n');
fprintf('DSB-SC: %.4f\n', corr_dsb(1,2));
fprintf('SSB-SC USB: %.4f\n', corr_ssb_usb(1,2));
fprintf('SSB-SC LSB: %.4f\n', corr_ssb_lsb(1,2));

%% Additional Figure: Detailed Spectrum Analysis
figure('Position', [150, 150, 1000, 600]);

subplot(2,2,1);
plot(f/1000, 20*log10(abs(M_f)+eps), 'b', 'LineWidth', 1.5);
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
title('Message Signal Spectrum (dB)');
grid on;
xlim([-5 5]);

subplot(2,2,2);
plot(f/1000, 20*log10(abs(DSB_f)+eps), 'g', 'LineWidth', 1.5);
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
title('DSB-SC Spectrum (dB)');
grid on;
xlim([-fc/1000*1.5 fc/1000*1.5]);

subplot(2,2,3);
plot(f/1000, 20*log10(abs(SSB_USB_f)+eps), 'm', 'LineWidth', 1.5);
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
title('SSB-SC (USB) Spectrum (dB)');
grid on;
xlim([-fc/1000*1.5 fc/1000*1.5]);

subplot(2,2,4);
plot(f/1000, 20*log10(abs(SSB_LSB_f)+eps), 'c', 'LineWidth', 1.5);
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
title('SSB-SC (LSB) Spectrum (dB)');
grid on;
xlim([-fc/1000*1.5 fc/1000*1.5]);

sgtitle('Frequency Domain Analysis (dB Scale)', 'FontSize', 14, 'FontWeight', 'bold');
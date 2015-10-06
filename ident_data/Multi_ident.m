%% Data Preprocessing
%
clc,

cd('/home/lex/Documents/Matlab_Files/ident_data/Sawtooth Waveform')
input_swt = load('input_sawtooth.csv');
output_swt = load('output_sawtooth.csv');

cd('/home/lex/Documents/Matlab_Files/ident_data/MLS Sequence');
input_mls = load('input_MLS.csv');
output_mls = load('output_MLS.csv');

cd('/home/lex/Documents/Matlab_Files/ident_data/UWGN');
input_ugwn = load('inputUGWN.csv');
output_ugwn = load('outputUGWN.csv');

save sawtoothdata.mat output_swt input_swt;
save mlsdata.mat output_mls input_mls;
save ugwndata.mat output_ugwn input_ugwn;

load sawtoothdata.mat;
load mlsdata.mat;
load ugwndata.mat;

clf;
format compact

%% Compute the FFTs and Power Spectrum of each input signal.
fft_mls    = fft(input_mls);
fft_ugwn   = fft(input_ugwn);
fft_swt    = fft(input_swt);

% Compute the two-sided spectrum P2 then compute the single-sided P1 based
% on P2 and even-valued signal length L
L_mls = length(input_mls);
L_ugwn = length(input_ugwn);
L_swt = length(input_swt);

T = (167/3) / 1000;         %sampling period
Fs = 1/T;                   %sampling frequency

%mls compute
P2_mls = abs(fft_mls/L_mls);
P1_mls = P2_mls(1:L_mls/2+1);
P1_mls(2:end-1) = 2*P1_mls(2:end-1);

f_mls = Fs*(0:(L_mls/2))/L_mls;
figure(1), subplot(211)
plot(f_mls,P1_mls, 'LineWidth', 1.4)
title('Spectrum of MLS Signal. Poly Order = 5')
xlabel('f (Hz)')
ylabel('|P1_{mls}(f)|')

%ugwn compute
P2_ugwn = abs(fft_ugwn/L_ugwn);
P1_ugwn = P2_ugwn(1:L_ugwn/2+1);
P1_ugwn(2:end-1) = 2*P1_ugwn(2:end-1);

f_ugwn = Fs*(0:(L_ugwn/2))/L_ugwn;
subplot(212)
plot(f_ugwn,P1_ugwn, 'LineWidth', 1.4)
title('Spectrum of Uniform White Gaussian Noise Signal')
legend('\sigma = 165/\sqrt(1.732)')
xlabel('f (Hz)')
ylabel('|P1_{ugwn}(f)|')

figure(2), subplot(211)
plot(L_mls, input_mls, 'LineWidth', 1.4)
title('Time series plot of PRBS Maximum Length Sequence Signal, ')
xlabel('time (msecs)'), ylabel ('current (mA) ')

subplot(212)
plot(L_ugwn, input_ugwn, 'LineWidth', 1.4)
title('Time series plot of Uniform Gaussian White Noise Signal')
xlabel('time (msecs)'), ylabel ('current (mA) ')

%Sawtooth compute
P2_swt = abs(fft_swt/L_swt);
P1_swt = P2_swt(1:L_swt/2+1);
P1_swt(2:end-1) = 2*P1_swt(2:end-1);

f_swt = Fs*(0:(L_swt/2))/L_swt;
figure(3)
plot(f_swt,P1_swt)
title('Single-Sided Amplitude Spectrum of Sawtooth Input Signal')
xlabel('f (Hz)')
ylabel('|P1_{swt}(f)|')

%% Raw Time-Series Signal
figure(4), subplot(211)
plot(output_mls, 'LineWidth', 1.4)
title('Raw MLS TIme Series Data'),
legend('Sampling Frequency, Fs = 3 X Fusion Freq.')
xlabel('Time(mSecs)'), ylabel('Head Position (mm)')

subplot(212)
plot(input_mls, 'LineWidth', 1.4)
title('Raw MLS Time Series Data'),
legend('Sampling Frequency, Fs = 3 X Fusion Freq.'),
xlabel('Time(mSecs)'), ylabel('Current (mA)')

%% Detrending Input - Output Data
input_mlsd = detrend(input_mls);
output_mlsd = detrend(output_mls);
figure(5), subplot(211)
plot(output_mlsd, 'LineWidth', 1.4)
title('Detrended MLS TIme Series Data'),
legend('Sampling Frequency, Fs = 3 X Fusion Freq.')
xlabel('Time(mSecs)'), ylabel('Head Position (mm)')

subplot(212)
plot(input_mlsd, 'LineWidth', 1.4)
title('Detrended MLS Time Series Data'),
legend('Sampling Frequency, Fs = 3 X Fusion Freq.'),
xlabel('Time(mSecs)'), ylabel('Current (mA)')

save mlsdetrended.dat input_mlsd output_mlsd
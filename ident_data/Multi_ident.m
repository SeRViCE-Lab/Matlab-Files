%% Data Preprocessing
%
clc, clear
close all

format compact;

Ts = 1;              %Sampling time for all three expts

cd('/home/lex/Documents/Matlab_Files/ident_data/Sawtooth Waveform')
input_swt = load('input_sawtooth.csv');
output_swt = load('output_sawtooth.csv');

swt = iddata(input_swt, output_swt, Ts);

cd('/home/lex/Documents/Matlab_Files/ident_data/MLS Sequence');
input_mls = load('input_MLS.csv');
output_mls = load('output_MLS.csv');

mls = iddata(input_mls, output_mls, Ts);

cd('/home/lex/Documents/Matlab_Files/ident_data/UWGN');
input_ugwn = load('inputUGWN.csv');
output_ugwn = load('outputUGWN.csv');

ugwn = iddata(input_ugwn, output_ugwn(1:length(input_ugwn)), Ts);

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
P1_ugwn(2:end-1)= 2*P1_ugwn(2:end-1);

f_ugwn = Fs*(0:(L_ugwn))/L_ugwn;

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
%% UWGN Plots
clf, clc

subplot(211),
t_ugwn = T * linspace(1, L_ugwn, L_ugwn);
plot(t_ugwn, 0.5* input_ugwn, 'LineWidth', 0.6),
title('Uniform Gaussian White Noise Signal'),
legend('Sampling Frequency = 17.9641Hz'),
xlabel('Time (Seconds)'), ylabel ('Current (mA) ')

subplot(212),
plot(f_ugwn(1:8969),P2_ugwn, 'LineWidth', 1.8)
title('Spectrum of Uniform White Gaussian Noise Signal')
legend('\sigma = 165/sqrt(3)')
xlabel('f (Hz)'),  ylabel('Normalized UWGN FFT'),
ylim([-0.2, 5])
%% Raw MLS Time-Series  Signal
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

%% Raw UWGN Time-Series  Signal
clf, subplot(211)
plot(output_ugwn, 'LineWidth', 1.4)
title('KF Fusion Estimate of both Sensor''s Observation'),
legend('Sampling Frequency, Fs = 3 X Fusion Freq.')
xlabel('Time(mSecs)'), ylabel('Head Position (mm)')

subplot(212)
plot(input_ugwn, 'LineWidth', 0.8)
title('Input Excitation Signal'),
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
%% LQG test
clear, clc, close all
cd('/home/lex/Documents/Matlab_Files/ident_data/MLS Sequence')
load('FreeForm.mat');
ss_test1

% setup
A=ss_test1.a;
B=ss_test1.b;
C=ss_test1.c;


Q=[1000000, 0;
    0, 100];
    
R= 0.1;

[K,S,E] = dlqr(A,B,10000*Q, R),
[n, d]  = ss2tf(A-B*K,B,C, 0);

DC=sum(n)/sum(d), %dc gain found evaluating TF at z=1

%fprintf('        [Press ANY key to continue...]\n\n');
%pause,

%clf
sysol = ss(A,B,C,0, ss_test1.Ts);
subplot(211)
step(sysol)
title('open loop response', 'LineWidth', 2.5)
grid

%fprintf('        [Press ANY key to continue...]\n\n');
%pause,

%clf
syscl = ss(A-B*K,1/DC*B,C,0,1/11);
subplot(212)
step(syscl)
title('closed loop', 'LineWidth', 2.5)
grid

%% Oct 17 Expt 1.I Plots
close all; clc;
cd('/home/lex/Documents/Matlab_Files/Expt_Plots/Expt I')
input = load('current.csv');
yhat  = load('yhat.csv');
y     = load('fusion.csv');
ref   = load('setpoint.csv');

format compact
clf

time = 0.1667 * linspace(1, length(y), length(y));
subplot(221),
p = plot(time, [ref, y, yhat]),
hold on,

p(1).LineWidth = 2.5;
p(2).LineWidth = 2.5;
p(3).LineWidth = 2.5;

xlabel('Time (seconds)'), 
ylabel('Head Position Relative to Camera Center (Millimeters)'),

legend('Reference', 'Kinect Fused yhat(kT)', 'LQG Estimated yhat(kT)', 'location', 'best'),

title('Experiment I'),

xlim([15, 100]), ylim([684, 710])

hold off

% Expt 1.II Plots
 clc;
input1 = load('current1.csv');
yhat1  = load('yhat1.csv');
y1     = load('fusion1.csv');
ref1   = load('setpoint1.csv');

%clf

time1 = 0.1667 * linspace(1, length(y1), length(y1));
subplot(222),
p1 = plot(time1, [ref1, y1, yhat1]),
hold on,

p1(1).LineWidth = 2.5;
p1(2).LineWidth = 2.5;
p1(3).LineWidth = 2.5;

xlabel('Time (seconds)'), 
ylabel('Head Position Relative to Camera Center (Millimeters)'),

legend('Reference', 'Kinect Fused yhat(kT)', 'LQG Estimated yhat(kT)', 'location', 'best'),

title('Experiment II'),

xlim([10, 107]), ylim([680, 710])

hold off

% Expt 1.III Plots
 clc;
input2 = load('current2.csv');
yhat2  = load('yhat2.csv');
y2     = load('fusion2.csv');
ref2   = 693.47 * ones( length(y2), 1);

%clf

time2 = 0.1667 * linspace(1, length(y2), length(y2));
subplot(223)
p2 = plot(time2, [ref2, y2, yhat2]),
hold on,

p2(1).LineWidth = 2.5;
p2(2).LineWidth = 2.5;
p2(3).LineWidth = 2.5;

xlabel('Time (seconds)'), 
ylabel('Head Position Relative to Camera Center (Millimeters)'),

legend('Reference', 'Kinect Fused yhat(kT)', 'LQG Estimated yhat(kT)', 'location', 'best'),
title('Experiment III'),

xlim([10, 90]), ylim([680, 695]);

hold off

% Expt 1.IV Plots
 clc;
input3 = load('current3.csv');
yhat3 = load('yhat3.csv');
y3     = load('fusion3.csv');
ref3   = load('setpoint3.csv');

%clf

time3 = 0.1667 * linspace(1, length(y3), length(y3));
subplot(224),
p3 = plot(time3, [ref3, y3, yhat3]),
hold on,

p3(1).LineWidth = 2.5;
p3(2).LineWidth = 2.5;
p3(3).LineWidth = 2.5;

xlabel('Time (seconds)'), 
ylabel('Head Position Relative to Camera Center (Millimeters)'),

legend('Reference', 'Kinect Fused yhat(kT)', 'LQG Estimated yhat(kT)', 'location', 'best'),

title('Experiment IV'),

xlim([10, 100]), ylim([685, 700]);

hold off

%% Oct 18 Expts
% Expt I Plots
close all; clc;
cd('/home/lex/Documents/Matlab_Files/Expt_Plots/Expt II')
input = load('current.csv');
yhat  = load('LQGyhat.csv');
y     = load('Fusionyhat.csv');
ref   = load('setpoint.csv');

format compact
clf

time = 0.1667 * linspace(1, length(y), length(y));
subplot(221),
p = plot(time, [ref, y, yhat]),
hold on,

p(1).LineWidth = 2.5;
p(2).LineWidth = 2.5;
p(3).LineWidth = 2.5;

xlabel('Time (seconds)'), 
ylabel('Head Position Relative to Camera Center (Millimeters)'),

legend('Reference', 'Kinect Fused yhat(kT)', 'LQG Estimated yhat(kT)', 'location', 'best'),

title('Experiment I'),

xlim([15, 145]), ylim([684, 695])

hold off

% Expt II Plots
 clc;
input1 = load('current1.csv');
yhat1  = load('LQGyhat1.csv');
y1     = load('Fusionyhat1.csv');
ref1   = load('setpoint1.csv');

%clf

time1 = 0.1667 * linspace(1, length(y1), length(y1));
subplot(222),
p1 = plot(time1, [ref1, y1, yhat1]),
hold on,

p1(1).LineWidth = 2.5;
p1(2).LineWidth = 2.5;
p1(3).LineWidth = 2.5;

xlabel('Time (seconds)'), 
ylabel('Head Position Relative to Camera Center (Millimeters)'),

legend('Reference', 'Kinect Fused yhat(kT)', 'LQG Estimated yhat(kT)', 'location', 'best'),

title('Experiment II'),

xlim([10, 98]), ylim([680, 700])

hold off

% Expt III Plots
 clc;
input2 = load('current2.csv');
yhat2  = load('LQGyhat2.csv');
y2     = load('Fusionyhat2.csv');
ref2   = load('setpoint2.csv');

%clf

time2 = 0.1667 * linspace(1, length(y2), length(y2));
subplot(223)
p2 = plot(time2, [ref2, y2, yhat2]),
hold on,

p2(1).LineWidth = 2.5;
p2(2).LineWidth = 2.5;
p2(3).LineWidth = 2.5;

xlabel('Time (seconds)'), 
ylabel('Head Position Relative to Camera Center (Millimeters)'),

legend('Reference', 'Kinect Fused yhat(kT)', 'LQG Estimated yhat(kT)', 'location', 'best'),
title('Experiment III'),

xlim([10, 98]), ylim([680, 695]);

hold off

% Expt IV Plots
 clc;
input3 = load('current3.csv');
yhat3 = load('LQGyhat3.csv');
y3     = load('Fusionyhat3.csv');
ref3   = load('setpoint3.csv');

%clf

time3 = 0.1667 * linspace(1, length(y3), length(y3));
subplot(224),
p3 = plot(time3, [ref3, y3, yhat3]),
hold on,

p3(1).LineWidth = 2.5;
p3(2).LineWidth = 2.5;
p3(3).LineWidth = 2.5;

xlabel('Time (seconds)'), 
ylabel('Head Position Relative to Camera Center (Millimeters)'),

legend('Reference', 'Kinect Fused yhat(kT)', 'LQG Estimated yhat(kT)', 'location', 'best'),

title('Experiment IV'),

xlim([20, 100]), ylim([680, 695]);

hold off
%% Feedforward Gain Computation
clear all; clc
cd('/home/lex/Documents/Matlab_Files/ident_data/UWGN')

uwgn_ident = load('uwgn_ident.mat');  % load the identification model
uwgn_ident = uwgn_ident.uwgn;

% Construct state space from ident data
A = uwgn_ident.a; B = uwgn_ident.b; C = uwgn_ident.c; D = uwgn_ident.d;
Ts = uwgn_ident.Ts;

Sys = ss(A, B, C, D, Ts);

%LQR Model
Qy = diag(1.51991);
Ry = diag(0.862916);
Ny = [0]';

Q = diag([1.51991, 8.11131]);
R = diag(0.862916);
N = [0, 0]';

[Ky, Sy, ey] = lqry(Sys, Qy, Ry, Ny);
[K, S, e] = lqr(A, B, Q, R, N);

%feedforward gain $k_g$
Kt = [-2.50099, -2.14999];
inner = -pinv(A - B * Kt);
k_g = pinv(C * inner * B)


%Retrieve the dcgain of the model and fetch feedforward gain
%Kff = 1/dcgain(uwgn)

%% Design of Model Algorithmic Controller
clear all; clc
cd('/home/lex/Documents/Matlab_Files/ident_data/UWGN')

uwgn_data = load('uwgn_data.mat');
uwgn_data = uwgn_data.uwgn_data;
uwgn_ident = load('uwgn_ident.mat');  % load the identification model
uwgn_ident = uwgn_ident.uwgn;

% We convert the obtained model to FIR polynomial model eqn 2.2 Soeterboek
uwgn_poly = idpoly(uwgn_ident);

%{ 
% this is unneeded as covariance model is preserved in polynomial sampling
%use zero iteration to recompute lost covariance
opt = polyestOptions;
opt.SearchOption.MaxIter = 0;
Ts = 0.0556;

%convert input-output time-domain signals into iddata object
uwgn_iddata = iddata(uwgn_data.output_ugwn, uwgn_data.input_ugwn, Ts);
uwgn_ARMAX = polyest(uwgn_iddata, uwgn_poly, opt);

fcn = @(x)idpoly(x);
uwgn_ARMAXt = translatecov(fcn, uwgn_ARMAX);
%}


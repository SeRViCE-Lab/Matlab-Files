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
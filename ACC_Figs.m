clc, close all
cd('/home/lex/Documents/Matlab_Files/kalman_data/kman_data6')
obskf1 = load('Xbox_obs2.csv') - 180;
predkf1 = load('Xbox_Pred2.csv') - 180;
updatekf1 = load('Xbox_Updates2.csv') - 180;
    
subplot(211)
plot(obskf1, 'LineWidth', 1.7)
xlabel('Time(sec)'); ylabel('Depth (mm)'), grid on
title('Kinect Xbox Tracking Noise on a Static Target at 684mm')
xlim([0, 2300])
legend('Ground Truth = 680mm; Covariance = 21.5022 ',  'location', 'northeast')
hold off

%Protonect
cd('/home/lex/Documents/Matlab_Files/')
points = load('depthPoints.csv');
time = linspace(1, length(points), length(points));
subplot(212)
plot(time, points, 'LineWidth', 1.5); hold on;
xlabel('Time(sec)'); ylabel('Depth (mm)');
ylim([994 1010]); xlim([100 1000]);
legend('Ground Truth = 960mm; Covariance = 3.1912', 'location', 'northeast')
title('Kinect v1  Position Estimation of a Static Object');
hold off;
%%
cd('/home/lex/Documents/Matlab_Files/');
rawdepth = load('ROSFacePoints.csv');
pred = load('ROSPrediction.csv');
corr = load('ROSCorrected.csv');
update = load('ROSUpdates.csv');
subplot (212)
plot( linspace(1, length(rawdepth), length(rawdepth)), rawdepth, 'LineWidth', 1.5),
hold on
xlabel('Samples'), ylabel('Depth (mm)'),
%ylim([670 690]), xlim([0 192]);
title('Kinect v1 Position Estimation of a static object')
legend('Ground Truth = 680mm', 'location', 'northeast')
hold off

%% Expt 1 Protonect
figure(1)
cd('/home/lex/Documents/Matlab_Files/kalman_data')
measure = load('rawmeasures.csv');
time = linspace(1, length(measure), length(measure));
subplot(211)
plot(time, measure, 'LineWidth', 1.5); hold on;
xlabel('Time(sec)'); ylabel('Depth (mm)');
legend('Ground Truth = 680mm; Covariance = 22.7057 ',  'location', 'northeast')
xlim([0 340]); ylim([665 695]);
title('Kinect X-box  Position Estimation of a static object.');
hold off;

%Protonect
cd('/home/lex/Documents/Matlab_Files/')
points = load('depthPoints.csv');
time = linspace(1, length(points), length(points));
subplot(212)
plot(time, points, 'LineWidth', 1.5); hold on;

%plot(time, eye_center, 'LineWidth', 2.5),
xlabel('Time(sec)'); ylabel('Depth (mm)');
ylim([994 1010]); xlim([0 1050]);
legend('Ground Truth = 960mm; Covariance = 3.1912', 'location', 'northeast')
title('Kinect v1  Position Estimation of a Static Object');
hold off;

%% ROS DEPTH?KALMAN
cd('/home/lex/Documents/Matlab_Files/');
rawdepth = load('ROSFacePoints.csv');
pred = load('ROSPrediction.csv');
corr = load('ROSCorrected.csv');
update = load('ROSUpdates.csv');
figure(2)
plot( linspace(1, length(rawdepth), length(rawdepth)), rawdepth, 'LineWidth', 1.5),
hold on
xlabel('Samples'), ylabel('Depth (mm)'),
%ylim([670 690]), xlim([0 192]);
title('Kinect v1 Position Estimation of a static object')
legend('Ground Truth = 680mm', 'location', 'northeast')
hold off
%% Kalman Results
close all
cd('/home/lex/Documents/Matlab_Files/kalman_data/kman_data3')
rawdepth = load('ROSFacePoints.csv');
pred = load('ROSPrediction.csv');
update = load('ROSUpdates.csv');

subplot(211)
plot( linspace(1, length(rawdepth), length(rawdepth)), rawdepth, 'LineWidth', 2.5),
hold on
xlabel('Time(Secs'), ylabel('Depth in mm'), grid on
xlim([1, 450])
title('Kinect v1 Tracking Noise on a Static Target')
legend('Observation', 'location', 'best')
hold off

subplot(212)
p = plot(linspace(1, length(pred), length(pred)), [pred, update(1:length(pred))]),
p(1).LineWidth = 1.5;
p(1).Marker    = '*';
p(1).MarkerFaceColor = 'r';

p(2).LineWidth = 6.5;
p(2).Marker    = '.';
p(2).MarkerFaceColor = 'b';
%p(2).MarkerEdgeColor = 'b';

grid on

hold on
xlabel('Time(secs)'), ylabel('Kalman Filter Estimates(mm)')
legend('Priori', 'Posteriori', 'location', 'best')
xlim([1, 450])
title('Recursive KF Results of Raw Observation in Real Time'), 
%% Adjusted covariance for Xbox 52.6832
clc, close all
cd('/home/lex/Documents/Matlab_Files/kalman_data/kman_data6')
obskf1 = load('Xbox_obs2.csv') - 180;
predkf1 = load('Xbox_Pred2.csv') - 180;
updatekf1 = load('Xbox_Updates2.csv') - 180;
    
subplot(211)
plot(obskf1, 'LineWidth', 1.7)
xlabel('Samples'), ylabel('Observation in mm'), grid on
title('Kinect Xbox Tracking Noise on a Static Target at 684mm')
xlim([100, 2300])
legend('Observation', 'location', 'northeast')
hold off

subplot(212)
p(2) = plot( linspace(1, length(updatekf1), length(updatekf1)), updatekf1, 'LineWidth', 1.5);
%p(1).LineWidth = 2.5;
%p(1).MarkerFaceColor = 'y';

p(2).LineWidth = 4.5;
p(2).MarkerEdgeColor = 'b';
legend('Estimates',...
     'location', 'best'),
xlim([100, 2300])
grid on, hold on
xlabel('Samples'), ylabel('Kalman updates in mm')
title('Kalman Filter Estimates with Kinect Xbox. [Truth: 684mm, R: 70mm^2]'),
hold off

%% SNR for both local KFs
cd('/home/lex/Documents/Matlab_Files/kalman_data/kman_data6')
obskf1 = load('Xbox_obs2.csv') - 180;
predkf1 = load('Xbox_Pred2.csv') - 180;
updatekf1 = load('Xbox_Updates2.csv') - 180;
snr_xbox = snr(updatekf1, obskf1)

cd('/home/lex/Documents/Matlab_Files/');
rawdepth = load('ROSFacePoints.csv');
pred = load('ROSPrediction.csv');
corr = load('ROSCorrected.csv');
update = load('ROSUpdates.csv');
snr_kv1 = snr(update, rawdepth(1:191))

%% Variance-weighted Fusion 1
clc; close all; clear all;
cd('/home/lex/Documents/Matlab_Files/kalman_data/kman_data_fused')
xkk1= load('Xbox_update.csv') - 200;
xkk2 = load('ros_estimates.csv') - 200;
xkk = load('easy_fusion.csv') - 200;
    
figure();
subplot(211)
plot( linspace(1, length(xkk1), length(xkk1)), [xkk1, xkk2, xkk], 'LineWidth', 4.5),
legend('Xbox Estimates', 'Kinect2 Estimates', 'Variance-Weighted Estimates',...
     'location', 'best'),
grid on, hold on
xlabel('Seconds'), ylabel('KF updates(mm)')
xlim([5 65])
title('Fusion of KF Estimate Tracks from Kinects Xbox and v1 . [Truth = 655 mm, Q = 1500mm^2, R_{Xbox} = 70mm^2, R_{K_2} = 70mm^2]'),
hold off

% Variance-weighted Fusion 2
clc; 
cd('/home/lex/Documents/Matlab_Files/kalman_data/kman_data_fused/kman_data_fused2')
xkk11= load('Xbox_estimates.csv');
xkk21 = load('ros_estimates.csv');
xkk1 = load('easy_fusion.csv');
    
subplot(212)
plot( linspace(1, length(xkk11), length(xkk11)), [xkk11, xkk21, xkk1], 'LineWidth', 3.5),
legend('Xbox Estimates', 'Kinect2 Estimates', 'Variance-Weighted Estimates',...
     'location', 'best'),
grid on, hold on
xlabel('Seconds'), ylabel('KF updates (mm)')
xlim([0 400]), %ylim([670 710])
title('Fusion of KF Estimate Tracks from Kinects Xbox and v1. [Truth = 687 mm, Q = 1500mm^2, R_{Xbox} = 100mm^2, R_{K_2} = 70mm^2]')
hold off
%% Variance-weighted Fusion 3; Q now 2000, R1 = 100, R2 = 60
clc; close all
cd('/home/lex/Documents/Matlab_Files/kalman_data/kman_data_fused/kman_data_fused2')
xkk12= load('Xbox_estimates2.csv');
xkk22 = load('ros_estimates2.csv');
xkk222 = load('easy_fusion2.csv');
    
figure();
plot( linspace(1, length(xkk12), length(xkk12)), [xkk12, xkk22, xkk222], 'LineWidth', 4.5),
legend('Xbox Estimates', 'Kinect2 Estimates', 'Variance-Weighted Estimates',...
     'location', 'best'),
grid on, hold on
xlabel('Seconds'), ylabel('Kalman estimates in mm')
title('Fusion of KF Estimate Tracks from Kinects Xbox and v1. [Truth = 687 mm, Q = 2000mm^{2}, Rxbox = 100, RK2 = 70mm^{2}]'),
hold off

%figure()
%plot(xkk22)



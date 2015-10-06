clc, close all
cd('/home/lex/Documents/Matlab_Files/')
points = load('depthPoints.csv');
time = linspace(1, length(points), length(points));
eye_center = 997 * ones(length(time));
plot(time, points, 'LineWidth', 1.5); hold on;

%plot(time, eye_center, 'LineWidth', 2.5),
xlabel('Time(sec)'); ylabel('Depth measurements (mm)');
xlim([0 900]);
title('Static Head/Right Eye center at 960mm from camera');
hold off;

%% Experiment 4: Kalman Filter and Baye's Measurement Updating
clc; figure, close all;
cd('/home/lex/Documents/Matlab_Files/');

%load data
data_one = load('kalmanpoints_1.csv');
data_two = load('kalmanpoints_2.csv');

time4 = linspace(1, length(data_one), length(data_one));
%{
coeffMA_one = ones(1, floor(max(data_one)))...
                    /floor(max(data_one)) ;
coeffMA_oner = 1000*coeffMA_one;
avgVal = filter(coeffMA_one, 1, data_one);
%}

cubicMA = sgolayfilt(data_one, 3, 201);
quarticMA = sgolayfilt(data_one, 4, 201);
quinticMA = sgolayfilt(data_one, 5, 901);

%plot(time4, data_one); hold on;
figure;
plot(time4, [data_one cubicMA quarticMA quinticMA ]); hold on;
%plot(time4,  avgVal);
xlabel('Time(msec)'); ylabel('Depth measurements (mm)');
xlim([0, 850]);
title('Depth measurements on a tracked object at 960(mm)');
%legend('Kinect Measurements', 'location', 'best');
legend('Kinect Measurements', 'Cubic-Weighted MA', 'Quartic-Weighted MA', ...
       'Quintic-Weighted MA','location', 'best');

hold off;

%% Friday Idea
close all; 
cd('/home/lex/Documents/Matlab_Files/');
data_one = load('kalmanpoints_1.csv');
data_frid = data_one(1:200); 
cubicMATest = sgolayfilt(data_frid, 4, 199); 
[m, n] = size(data_one);
    plot(data_frid), hold on
    plot (cubicMATest, 'LineWidth', 3.5);
    title('Kinect ToF Camera (Protonect) , Stationary object at 1002mm'); grid;
    xlabel('Samples'), ylabel('Distance from cam center (mm)');
    linkdata on;
    legend('Measurements', 'Savitzky-Golay Smoothed', 'location', 'best') 
for k = 201:m
    data_frid(k,:) = data_one(k);
    cubicMATest    = sgolayfilt(data_frid, 4, 199);
end
%pause(0.5)   
hold off; pause off

%% August 06
clc; close all;
cd('/home/lex/Documents/Matlab_Files/kalman_data')
measure = load('rawmeasures.csv');
time = linspace(1, length(measure), length(measure));
facepoints = 680 * ones(length(time));
subplot(311)
plot(time, measure, 'b-.'); hold on;
xlabel('Time(sec)'); ylabel('Depth measurements (mm)');
legend('Raw Measurements', 'location', 'best')
%ylim([980 1100]), xlim([0 800]);
title('Raw Measurements from Kinect X-box. Ground Truth = 680mm');
hold off;
%plot prediction values
predict = load('Prediction.csv');
time2 = linspace(1, length(predict), length(predict));
subplot(312)
plot(time2, [predict], 'r-.'); hold on;
xlabel('Time(sec)'); ylabel('Depth measurements (mm)');
legend('Kalman Prediction', 'True Head Position')
%ylim([980 1100]), xlim([0 800]);
title('Kalman Prediction with Kinect X-box');
hold off;
%plot updates
updates = load('Updates.csv');
time3 = linspace(1, length(updates), length(updates));
subplot(313)
plot(time3, [updates], 'g-.'); hold on;
xlabel('Time(sec)'); ylabel('Depth measurements (mm)');
legend('Kalman Updates', 'True Head Position')
ylim([670 700]), xlim([0 350]);
title('Kalman Updates with Kinect X-box');
hold off;
%% ROS DEPTH?KALMAN
close all
cd('/home/lex/Documents/Matlab_Files/');
rawdepth = load('ROSFacePoints.csv');
pred = load('ROSPrediction.csv');
corr = load('ROSCorrected.csv');
update = load('ROSUpdates.csv');
%{
%check for out of ordinary values
for a=1:length(rawdepth)
    if rawdepth(a) > 690
        rawdepth(a) = rawdepth(a-1)
    end
end

for b=1:length(corr)
    if corr(b) > 726
        corr(b) = corr(b-1)
    end
end
%}

figure()
subplot(211)
plot( linspace(1, length(rawdepth), length(rawdepth)), rawdepth, 'LineWidth', 1.5),
hold on
xlabel('Samples'), ylabel('Depth in mm')
title('raw depth values')
legend('With Kinect2', 'location', 'best')
hold off

subplot(212)
plot(linspace(1, length(pred), length(pred)), pred, 'LineWidth', 1.5),
hold on
xlabel('Samples'), ylabel('Kalman prediction in mm')
title('kalman prediction'), legend('With Kinect2', 'location', 'best')
hold off

figure(2)
subplot(211)
plot(linspace(1, length(corr), length(corr)), corr, 'LineWidth', 1.5),
hold on
xlabel('Samples'), ylabel('Corrected measurements in mm')
legend('With Kinect2', 'location', 'best')
title('Corrected measurement values')

subplot(2,1,2)
plot( linspace(1, length(update), length(update)), update,'LineWidth', 1.5),
hold on
xlabel('Samples'), ylabel('Kalman updates in mm')
legend('With Kinect2', 'location', 'best'),
title('kalman updates'),
hold off
%% Monday August 17 Data kman_data3

close all
cd('/home/lex/Documents/Matlab_Files/kalman_data/kman_data3')
rawdepth = load('ROSFacePoints.csv');
pred = load('ROSPrediction.csv');
update = load('ROSUpdates.csv');

figure()
subplot(211)
plot( linspace(1, length(rawdepth), length(rawdepth)), rawdepth, 'LineWidth', 1.5),
hold on
xlabel('Samples'), ylabel('Depth in mm'),
xlim([1, 450])
title('raw depth values')
legend('With Kinect2', 'location', 'best')
hold off

subplot(212)
plot(linspace(1, length(pred), length(pred)), pred, 'LineWidth', 2.5),
hold on
xlabel('Samples'), ylabel('Kalman Filter Results in mm')
xlim([1, 450])
title('kalman prediction'), legend('With Kinect2', 'location', 'best')
hold off

figure(2)
subplot(211)
plot( linspace(1, length(rawdepth), length(rawdepth)), rawdepth, 'LineWidth', 1.5),
hold on
xlabel('Samples'), ylabel('Corrected measurements in mm')
legend('With Kinect2', 'location', 'best')
title('Corrected measurement values')

subplot(2,1,2)
plot( linspace(1, length(update), length(update)), update,'LineWidth', 2.5),
hold on
xlabel('Samples'), ylabel('Kalman updates in mm')
legend('With Kinect2', 'location', 'best'),
title('kalman updates'),
hold off

%% Estimates of an estimate
close all, clc;
cd('/home/lex/Documents/Matlab_Files/kalman_data/kman_data4')
%kalman filter 1
obskf1 = load('KF1Corrected.csv');
predkf1 = load('KF1Prediction.csv');
updatekf1 = load('KF1Updates.csv');

obskf2 = load('KF2Corrected.csv');
predkf2 = load('KF2Prediction.csv');
updatekf2 = load('KF2Updates.csv');

figure(1);
plot( linspace(1, length(updatekf1), length(updatekf1)), [updatekf1, updatekf2], 'LineWidth', 1.5),
legend('Kalman filter 1 estimates', 'Kalman filter 2 estimates',...
     'location', 'best'),
grid on, hold on
xlabel('Samples'), ylabel('Kalman updates in mm')
title('Kalman updates'),
hold off

figure(2);
plot(linspace(1, length(updatekf1), length(updatekf1)), [predkf1, predkf2, obskf2], 'LineWidth', 2.5),
legend('Prediction of KF1', 'Prediction of KF2', 'kalman 1 updates', 'location', 'best')
grid on, hold on
xlabel('Samples'), ylabel('Kalman Predictions in mm')
title('Kalman Predictions'),
hold off
%% Online Measurement covariance using Xbox
close all, clear all, clc;
cd('/home/lex/Documents/Matlab_Files/kalman_data/kman_data5')
%kalman filter 1
obskf1 = load('obs.csv');
predkf1 = load('pred.csv');
updatekf1 = load('est.csv');
    
figure(1);
plot( linspace(1, length(updatekf1), length(updatekf1)), [predkf1, updatekf1], 'LineWidth', 1.5),
legend('Prediction', 'Estimates',...
     'location', 'best'),
grid on, hold on
xlabel('Samples'), ylabel('Kalman updates in mm')
title('Kalman updates with online R (k) estimation. [Truth: 1059]'),
hold off
%% Determine xbox estimate covariance
close all, clear all, clc;
cd('/home/lex/Documents/Matlab_Files/kalman_data/kman_data6')
%kalman filter 1
obskf1 = load('Xbox_obs.csv');
predkf1 = load('Xbox_Pred.csv');
updatekf1 = load('Xbox_Updates.csv');
    
figure();
plot( linspace(1, length(updatekf1), length(updatekf1)), [predkf1, updatekf1], 'LineWidth', 1.5),
legend('Prediction', 'Estimates',...
     'location', 'best'),
grid on, hold on
xlabel('Samples'), ylabel('Kalman updates in mm')
title('Kalman updates with online R (k) estimation. [Truth: 1059]'),
hold off
%% Adjusted covariance for Xbox 52.6832
clc, close all
obskf1 = load('Xbox_obs2.csv') - 180;
predkf1 = load('Xbox_Pred2.csv') - 180;
updatekf1 = load('Xbox_Updates2.csv') - 180;
    
subplot(211)
plot(obskf1, 'LineWidth', 1.7)
xlabel('Time(Secs'), ylabel('Depth in mm'), grid on
title('Kinect Xbox Tracking Noise on a Static Target at 684mm')
legend('Observation', 'location', 'northeast')
hold off

subplot(212)
p = plot( linspace(1, length(updatekf1), length(updatekf1)), [predkf1, updatekf1], 'LineWidth', 1.5),
legend('Prediction', 'Estimates',...
     'location', 'best'),
grid on, hold on
xlabel('Samples'), ylabel('Kalman updates in mm')
title('Kalman updates with online R (k) estimation. [Truth: 856mm, Cov: 52.6832]'),
hold off

%% Adjusted covariance for Xbox 70
close all
obskf1 = load('Xbox_obs3.csv');
predkf1 = load('Xbox_Pred3.csv');
updatekf1 = load('Xbox_Updates3.csv');
    
subplot(211)
plot(obskf1)

subplot(212);
plot( linspace(1, length(updatekf1), length(updatekf1)), [predkf1, updatekf1], 'LineWidth', 1.5),
legend('Prediction', 'Estimates',...
     'location', 'best'),
grid on, hold on
xlabel('Samples'), ylabel('Kalman updates in mm')
title('Kalman updates with online R (k) estimation. [Truth: 856mm, Cov: 70]'),
hold off


%% Variance-weighted Fusion 
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
xlabel('Samples'), ylabel('Kalman estimates in mm')
title('Sensor fusion of KF Estimate Tracks from Kinects one and two . [Truth = 655 mm, R_{Xbox} = 70mm^2, R_{K_2} = 70mm^2]'),
hold off
%% Variance-weighted Fusion 
clc; close all; clear all;
cd('/home/lex/Documents/Matlab_Files/kalman_data/kman_data_fused')
xkk1= load('Xbox_update.csv') -200;
xkk2 = load('ros_estimates.csv')-200;
xkk = load('easy_fusion.csv')-200;
    
figure();
subplot(211);
plot( linspace(1, length(xkk1), length(xkk1)), [xkk1, xkk2, xkk], 'LineWidth', 1.5),
legend('Xbox Estimates', 'Kinect2 Estimates', 'Variance-Weighted Estimates',...
     'location', 'best'),
grid on, hold on
xlabel('Samples'), ylabel('Kalman estimates in mm')
title('Sensor fusion of Observations from Kinects one and two . [Truth = 850 mm]'),
hold off
%% Variance-weighted Fusion 2
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
xlabel('Samples'), ylabel('Kalman estimates in mm')
title('Sensor fusion of Observations from Kinects one and two . [Truth = 687 mm (27'')]'),
hold off
%% Variance-weighted Fusion 3; Q now 2000, R1 = 100, R2 = 60
clc; 
cd('/home/lex/Documents/Matlab_Files/kalman_data/kman_data_fused/kman_data_fused2')
xkk12= load('Xbox_estimates2.csv');
xkk22 = load('ros_estimates2.csv');
xkk = load('easy_fusion2.csv');
    
figure();
plot( linspace(1, length(xkk1), length(xkk1)), [xkk1, xkk22, xkk2], 'LineWidth', 1.5),
legend('Xbox Estimates', 'Kinect2 Estimates', 'Variance-Weighted Estimates',...
     'location', 'best'),
grid on, hold on
xlabel('Samples'), ylabel('Kalman estimates in mm')
title('Sensor fusion of Observations from Kinects one and two . [Truth = 687 mm (27''), Q = 2000, Rxbox = 100, RK2 = 60]'),
hold off



%%Plot Depth Points
% Load the dang points
clc, 
cd('/home/lex/Documents/Matlab_Files/')
points = load('depthPoints.csv');
time = linspace(1, length(points), length(points));
eye_center = 997 * ones(length(time));
subplot(3, 1, 1), grid
plot(time, points); hold on;

plot(time, eye_center, 'LineWidth', 2.5),
xlabel('Time(msec)'); ylabel('Depth measurements (mm)');
%ylim([980 1010]), xlim([0 800]);
title('Static Head/Right Eye center at 960mm from camera');
hold off;

%% Experiment 2
clc
cd('/home/lex/Documents/Matlab_Files/')
points2 = load('depthPoints2.csv');
time2 = linspace(1, length(points2), length(points2));
eye_center2 = 997 * ones(length(time2));
%subplot(3,1,2);
figure

hold on;
xlabel('Time(msec)'); ylabel('Depth measurements (mm)');
%ylim([980 1010]), xlim([0 720]);
title('Static Head/Right Eye center at 960mm from camera');
%plot(time2, eye_center2, 'LineWidth', 2.5); grid;
hold off;
%% Experiment 3
clc; clf
cd('/home/lex/Documents/Matlab_Files/')
points3 = load('depthPoints3.csv');
time3 = linspace(1, length(points3), length(points3));
eye_center3 = 997 * ones(length(time3));
subplot (3, 1, 3);
plot(time3, points3); hold on;
%xlabel('Time(msec)'); ylabel('Depth measurements (mm)');
legend('True Measurements', 'Compensated Ground Truth')
ylim([980 1100]), xlim([0 800]);
title('Static Head/Right Eye center at 960mm from camera');
plot(time3, eye_center3, 'LineWidth', 2.5);
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
data_one = load('kalmanpoints_1.csv');
data_frid = data_one(1:200); 
cubicMATest = sgolayfilt(data_frid, 4, 199); 
plot(data_frid), hold on
plot (cubicMATest, 'LineWidth', 3.5);
xlabel('Time(msec)'); ylabel('Depth measurements (mm)');
linkdata on;
[m, n] = size(data_one(201:end));
legend('Measurements', 'Savitzky-Golay Smoothed', 'location', 'best') ;
for k = 201:m
    data_frid(k,:) = data_one(k);
    cubicMATest    = sgolayfilt(data_frid, 4, 199);
end
hold off


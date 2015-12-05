%% Nov 14 Plots

close all; clc;
cd('/home/lex/Documents/Matlab_Files/Vicon_Data')
trans = load('midface.csv');

%We find the intersection of the four points(= center of face)
x = trans(:,1);
y = trans(:,2);
z = trans(:,3);

%Translations follow the right hand rule from head in mid room of high bay

format compact
clf

time = 0.1667 * linspace(1, length(trans), length(trans));

p = plot(time, [x, y, z]),
hold on,

p(1).LineWidth = 2.5;
p(2).LineWidth = 2.5;
p(3).LineWidth = 2.5;

xlabel('Time (seconds)'), 
ylabel('Head Translation (Millimeters)'),

legend('X-Translation', 'Y-Translation', 'Z-Translation', 'location', 'best'),

title('Static Head Position Captured with Vicon'),

%% Produce Indiv Plots
clf
subplot(311);

p1 = plot(time, x),  hold on,
p1.LineWidth = 2.5;

xlabel('Time (seconds)'), ylabel('Head Translation (Millimeters)'),
axis([700 2000 -115 -114 ])
legend('X-Translation', 'location', 'best'),

title('Static Head Position Captured with Vicon'),

%Y-Translation
subplot(312);

p2 = plot(time, y, 'r'),  hold on,
p2.LineWidth = 2.5;

xlabel('Time (seconds)'), ylabel('Head Translation (Millimeters)'),
axis([700 2000 68 70 ])
legend('Y-Translation', 'location', 'best'),

title('Static Head Position Captured with Vicon'),

%Z-Translation
subplot(313);

p2 = plot(time, z, 'm'),  hold on,
p2.LineWidth = 2.5;

xlabel('Time (seconds)'), ylabel('Head Translation (Millimeters)'),
axis([700 2000 248 249 ])
legend('Z-Translation', 'location', 'best'),

title('Static Head Position Captured with Vicon'),

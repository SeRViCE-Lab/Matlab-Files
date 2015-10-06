%% Rot Parser
cd('/home/lex/Downloads/Fei Tao/')
clc; clear all; close all
Rot  = load('rotPP.csv');

%% Parse rotation matrices
R = zeros(3,3);  %initialize
i = 1;
R_1 = Rot(1:3, 1:3);
 for j = 1:3:length(Rot)-3
       Rs(:,:,i)=R_1'*Rot(j:j+2,:);
       rpy(i,:)=tr2rpy( Rs(:,:,i));
        i=i+1;
 end       
    


%%  Plot RPY
clf
p = plot(rpy)
p(1).LineWidth = 2.5;
p(2).LineWidth = 2.5;
p(3).LineWidth = 2.5;

legend('roll', 'pitch', 'yaw', 'location', 'northwest');
xlabel('Samples'), ylabel('Degrees?'),
title('RPY Head Motion With CMU Code');
    
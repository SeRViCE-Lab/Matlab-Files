cd('/home/lex/Documents/Matlab_Files/ident_data/Sawtooth Waveform')
load('ss_test.mat')

%% setup
A=ss_test.a;
B=ss_test.b;
C=ss_test.c;


    
%% - lqr
close all
Q=[100, 0, 0, 0;
    0, 10, 0, 0;
    0,  0, 10, 0;
    0, 0, 0, 10];
    
R= .03;

[K,S,E] = dlqr(A,B,100*Q,R)
[n, d]=ss2tf(A-B*K,B,C, 0)

DC=sum(n)/sum(d) %dc gain found evaluating TF at z=1

figure
sysol = ss(A,B,C,0,1/11);
step(sysol)
title('open loop')
grid

figure
syscl = ss(A-B*K,1/DC*B,C,0,1/11);
step(syscl)
title('closed loop')
grid



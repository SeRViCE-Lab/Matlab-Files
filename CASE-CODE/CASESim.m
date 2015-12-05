%%Matlab script for sys identification expt in my CASE 2015 Paper
clc; 
% I have provided the .csv files in this dir
%change as you see fit
cd('C:\Users\opo140030\Desktop\Current Data\Paper Data');  

%Load Modeling Data from RIO
currentdata             = load('currentdata.csv');      % Load input current data
kinectdata              = load('kinectdata.csv');        % Load Kinect (Output) Data 

% Define input and output data for sys ident preprocessing
outputdatatest          = kinectdata(1:8800);
inputdatatest           = currentdata(1:8800);

%Load Validation Data from RIO
currentdataval          = load('currentdataval2.csv');      % Load input current data
kinectdataval           = load('kinectdataval.csv');        % Load Kinect (Output) Data 

% Define input and output data for sys ident preprocessing
outputdataval           = kinectdataval(1:8800);
inputdataval            = currentdataval(1:8800);

srsdata                 = [outputdatatest, inputdatatest];
srsval                  = [outputdataval, inputdataval];

figure(1)
% Display the raw data;
%
fprintf('\n');
fprintf('Step 1: Display the raw data (output/input)\n\n');

clf reset,

subplot(2,1,1),     plot(outputdatatest);   hold on

title('OUTPUT: Raw Head Measurements')

xlabel('Samples'), ylabel ('Head Measurement (mm)')

axis([0 9000 630 650]); grid; hold off;

subplot(2,1,2),     plot(inputdatatest)

grid;   title('INPUT: Applied Current ')
    
xlabel('Samples'), ylabel ('Current (mA)')

fprintf('        [Press ANY key to continue...]\n\n');

pause,
%
% Remove data means to see structure in data as dataset dimension is low
%Geoff Hinton, NIPS 2006
fprintf('Step 2: Remove data means and then data trend\n\n');

fprintf('   2.1  Remove data means\n\n');

srs_detrended       =   detrend(srsdata);

figure(2);

subplot(2,1,1),     plot(srs_detrended(:,1)),   hold on
subplot(2,1,1),     grid,   title('Detrended Output')
xlabel('Samples'),  ylabel('Head Measurement (mm)'), 
axis([0 9000 -8 10]),    hold off;

subplot(2,1,2),     plot(srs_detrended(:,2))
subplot(2,1,2),     grid, title('Detrended Input')
xlabel('Samples'),  ylabel ('Current (mA)')
fprintf('        [Press ANY key to continue...]\n\n');
pause,
%
% Remove data trends; Just to ensure outliers and HF signals do not
% dominate model%

fprintf('   2.2  Detrended data\n\n');
srs_mean_removed        =    detrend(srsdata,0);
clf reset,
subplot(2,1,1),         plot(srs_mean_removed(:,1)),    hold on
subplot(2,1,1),         grid,       title('Mean Removed Output')
xlabel('Samples'),      ylabel ('Head Measurement (mm)'), 
axis([0 9000 -8 10]),   hold off;

subplot(2,1,2),         plot(srs_mean_removed(:,2)),
subplot(2,1,2),         title('Mean Removed Input'),    grid
pause,
%
% Display the cross correlation function between 
% the raw input and output data;
%
fprintf('Step 3: Correlation analysis for raw input/output\n\n');
%% Display the cross correlation function between 
% the raw input and output data;
%
fprintf('        [Press ANY key to continue...]\n\n');
pause,
%
% Display the auto-correlation function the detrended input 
% and the prewhitended input;
%
fprintf('Step 4: Correlation analysis for detrnded\n');
fprintf('        and pre-whitened input and output\n\n');
input_length            = length(srsdata);
z1                      = srs_detrended(:, 1);
z2                      = srs_detrended(:, 2);
vrawout                 = input_length * var( srsdata(:, 1) );
vrawin                  = input_length * var( srsdata(:, 2));
varz1                   = input_length * var( z1 );
varz2                   = input_length * var( z2 );
ccf_raw                 = crosscorr(srsdata(:, 1) - mean( srsdata(:, 1) ), ...
                            srsdata(:, 2) - mean( srsdata(:, 2) ));
ccf_detrended           = crosscorr(z1 - mean(z1), z2 - mean(z2) );
ccf_raw                 = ccf_raw/sqrt (vrawout * vrawin);
ccf_detrended           = ccf_detrended/sqrt(varz1*varz2);
clf reset,

plot(ccf_raw, '-')

title('CCF[Raw Input-Output]("-.") ;  CCF[Detrended Input-Output](".")'),

xlabel('lags'),     grid,   legend('Raw Input CCF', 'Detrended Input CCF',...
                             'location', 'northeast'),

%% Auto - correlation function of input
% Display the auto-correlation function of the detrended input; 
%
fprintf('   4.1  Auto-corr analysis of detrended input\n\n');
acf_in                  = autocorr(z2,z2);
subplot(111),

plot(acf_in,'-','LineWidth', 1.5),
                      
xlabel('lags'),

title('Auto-Correlation Function of Detrended Input')
grid
fprintf('        [Press ANY key to continue...]\n\n');
pause,

% Display the auto-correlation function of the pre-whitened input; 
%
%% Pre-whiten signals and filter since input signal was not Gaussian & White
 fprintf('   (a)  Select model order for input prewhitening:\n'); 
  m                     =   input('        m= ');
  th                    =   ar(z2, m);
  model_ar              =   ar(z2, m);
  
  fprintf('\n');
  fprintf('   The filter model for input prewhitening:\n\n');
  fprintf('   '),   model_ar
  fprintf('\n');
  fprintf('        [Press ANY key to continue...]\n\n');
  
  pause,
  
  fprintf('   (b)  Prewhitened output and input\n\n');
  
  b                     =   [0  th(3, :)];
  prewhitened_opt       =   z1 + filter(b, 1, z1);
  
  % Use the prediction error model for model recontruction
  e                     =   pe(model_ar, z2);   
  clf reset,
  subplot(211),         plot(prewhitened_opt),  hold on,
                        title('Pre-whitened Output'), 
                        xlabel('Samples'), grid
                        axis([0 9000 -5 10]), hold off
                        
  subplot(212),         plot(e), hold on,   
                        title('Pre-whitened Input'),  
                        xlabel('Samples'), grid
                        
  fprintf('             [Press ANY key to continue...]\n\n');
  
  pause,
  
  fprintf('   (c)  Input whiteness test\n\n');
  
  subplot(111),
                        [e,s]=resid_acse(z2,th);
                        
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
%  Check model residues
  [e,r]                 =   resid(model_ar, z2);
fprintf('        [Press ANY key to continue...]\n\n');
pause,
%
% Display the cross correlation function the pre-whitened
% input and output; 
%
  fprintf('   4.3  Cross correlation function between\n');
  fprintf('        the prewhitened input and output\n\n');
  ccf_opt_ipt           =   crosscorr_acse(prewhitened_opt,e,'coeff');
  x1                    =   ccf_opt_ipt(296:320,1);
  time                  =   [0 : 24];
  
  subplot(111),
                        plot(time, x1, '-', time, x1, 'o'),       grid
                        title('CCF Between the Prewhitened Input and Output'), 
                        grid
                        xlabel('lags')
                        
  fprintf('        [Press ANY key to continue...]\n\n');
  
  pause,
  
  clear x1
  
  fprintf('        Are you satisfied with the results?\n');
  fprintf('        Please type a number to confirm!\n\n');
  fprintf('        0--YES; Other Numbers--NO!\n\n');
  reply=input('        Your choice is [0/others]: ');
  fprintf('\n\n');
  
  if reply          ~=  0
     fprintf('        *****  Repeat Step 4.2! *****\n\n');
     reply           = 1;
  end
%% Cross Spectral Analysis 
%  
fprintf('Step 5: Cross spectral analysis for\n');
fprintf('        detrended input and output\n\n');
cross_in            = [z2];  %No Periodic extension for input;
cross_out           = [z1]; % No Periodic extension for output;
subplot(111),

spectrum_acse(cross_in,cross_out);
fprintf('        [Press ANY key to continue...]\n\n');
pause,
fprintf(' \n');
%
%% System Identification 
%  
fprintf('Step 6: System identification\n\n');
reply               =   1;
while reply         ~=  0
  fprintf('   6.1  Model order selection\n\n');
   
  no=input('        Model order for output, no=');
  ni=input('        Model order for input, ni=');
  nn=input('        Model order for noise, nn=');
  nd=input('        Time delay for input, nd=');
  fprintf('\n\n');
  
  order             =   [no ni nn nd];
  th                =    armax(srs_detrended,order);
  model_armax       =    armax(srs_detrended,order);
  
  fprintf('\n');
  fprintf('   The armax model for prewhitened input-output:\n\n')
  fprintf('   '), model_armax;
  fprintf('\n');
  fprintf('        [Press ANY key to continue...]\n\n');
  
  pause,
  
  e                 =       resid(srs_detrended,th);
  fprintf('   6.2  Model reponses\n\n');
%
% Predicted output
%
  fprintf('   (a)  One-step-ahead predicted outputs\n\n');
  
  predicted_opt     =       z1 - e;
  
  time              =   [1:1:length(srsdata)]';
  
  subplot(111),
                        plot(time, z1, '-', time,predicted_opt, '*')
                        title('Measurement("-"); One-Step-Ahead Predicted Output("*")'),
  
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
%
% Deterministic output
%
  fprintf('   (b)  Model predicted outputs\n\n');
  deterministic_opt = idsim_acse(z2,th);
  
  subplot(111),
  
  plot(time,z1, '-', time, deterministic_opt, '+'),
  title('Measurement("-"); Model Predicted Output("+")'),
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
%
  fprintf('   (c)  Predicted residuals and errors\n');
  de            =   z1 - deterministic_opt;
  time          =  [1:1:length(srsdata)]';
  clf reset,
  subplot(211),     
                    plot(e);
                    
  subplot(211), 
                    title('One-Step-Ahead Prediction Residuals'),
                    
  subplot(212),
                    plot(de)
                    
  subplot(212),
                    title('Model Prediction Errors'),
                    
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
  fprintf('   (d)  Impulse response\n\n');
  
  z3                =  [4.0, zeros(1, length(srsdata) - 1)]';
  impulse_response  =  idsim(z3, th);
  impy              =   [impulse_response(1:25,1)];
  time              =   [0:1:24]';
  
  subplot(111),
  plot(time, impy, '+', time, ccf_opt_ipt(296:320,1), '-', ...
                     time, ccf_opt_ipt(296:320,1), 'o'),
  title('Impulse Response "+" and CCF "-" Between Prewhitened Output and Input'),
  fprintf('        [Press ANY key to continue...]\n\n');
  
  pause,
  
  fprintf('   6.3  Model validity test\n\n');
  
  subplot(111),
  [e,s]             =   resid(srs_detrended,th);
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
  
%  fprintf('        Display the model validity test again!\n\n\n');
%  subplot(111),
%  [e,r]=resid(model_armax,gas_detrended);
%  pause,
  %
  fprintf('Step 7  Is the model OK?\n\n');
  fprintf('        Please type a number to confirm!\n\n');
  fprintf('        0--YES\n');
  fprintf('        Other Numbers--NO\n\n');
  
  reply             =   input('        Your choice is [0/others]: ');
  fprintf('\n\n');
  if reply          ~=  0
     fprintf('        ***** Repeat Step 6.1! *****\n\n');
     reply          =   1;
  end
end

close;
%% Transient Response Analysis

clc; scrsz          = get(0,'ScreenSize');
clc;                close all
s                   = tf('s');

%sstestingdata was the response from procest function in the system id toolbox
%Feel free to use matlab scripts for your process estimation
G                   = sstestingdata;

figure(1);

bode(G),        hold on,     grid on,   axis ON,    hold off;

set(0,'CurrentFigure',2)

subplot(2,1,1),         grid, 

step(G), 

Su = stepinfo(step(G));

display(Su);

title('Step Response of Identified Model with System Identification')
% Magnitude of PID constants. I used the pidtoolbox to find optimal pid
% gains
%Kp = 547.658; Kd = 68.947; Ki = 0.00561; wc = 6.3069; Gc = Kp+Kd*s+Ki/s;
%Real PID Constants
Kp                  =   547.658; 
Kd                  =   -68.947; 
Ki                  =   -0.00561; 
wc                  =   6.3069; 
Gc                  =   Kp+Kd*s+Ki/s;

newG                =   series(G, Gc); 
sys                 = feedback(newG,1);

subplot(2,1,2),     grid,   step(sys)
title('Closed Loop Unity Feedback Step Response of SRS system with Incorporated PID')

S = stepinfo(sys);
display(S)

%% Find bode plot of detrended srs data in system identification tool box and import to workspace
%%
% Detrend data by removing means/linear trends
%zd = detrend(z);

% Estimate frequency response and spectrum by spectral analysis with freq
% dependent resolution
%zd_fdr = spafdr(zd,[]);
zd_fdr = freq_resp_data;
bode(zd_fdr), hold on, grid ON, axis ON, 
title('Bode Frequency Response of Detrended Input and Output Data'), hold off % Plot the bode response

% Filter frequency response data with an ideal frequency domain filter
f1                  = 0.00232; 
f2                  = 0.644; 
f3                  = 0.1644; 
f4                  = 7.62;  % define filter ranges
zd_fdr_f            = idfilt( zd_fdr, [f1, f2; f3, f4], 'FilterOrder', 5, 'Causal');

% Model Identification
P2                  = idproc('P2DI', 'InputName', 'Current'); % Create a 2nd-order, deadtimee, integrating process idproc object
P2.Tp2.min          = 5;            % Set the lower bound for the pole location
P2_wn               = pem(zd_fdr_f,P2);   % Identify the model using the frequency

% Provide information about the model
present(P2_wn)
% Take model to state space
P2_ss               = tf2ss(P2_wn);

% Find cut-off frequency
omega_c_wn          = 1/P2_wn.Tp2.value;
fc_wn               = omega_c_wn/(2*pi)

% Percent error in cutoff frequency between model and measured data
fc_err              = abs(((fc_wn - fc_mw)/fc_mw)*100)

% Compare the time-domain response of this model with the measured data
compare(zd,P2_wn)
%% Residuals Analysis
s = tf('s');
srsmodel = (0.0006 * (s-1.7137)/((s+0.01)*(s+0.1028)))*exp(-2*s)  % Transfer Function

e = resid(srsdata, srsmodel, 'FR'), hold on, 
grid, title('Model Residuals Frequency Response Against Soft Robot System Frequency Response'), hold off;

%% Control Analysis
clc, clear all,
clf reset
s = tf('s');
Ge = (-0.0006 * (s-1.7137)/((s+0.01)*(s+0.1028)))        %*exp(-2*s) ; % Transfer Function

Gdel = exp(-2*s);

Pdel = pade(Gdel, 2),

% Plant and delay 
G = series(Pdel, Ge), 

% Solved PI Controller
close all;
clf reset
h1 = stepplot(G, ':b'),
p1 = getoptions(h1),
p1.Title.String = 'Plant Open-Loop Step Response';
p1.TimeUnits= 'seconds',
setoptions(h1, p1)
stepinfo(G)
 
% Primary Loop
fprintf('Press any key to continue'),
pause,
clf reset,


C =  pid(3.79, 0.0344, 0), 
Gpiol = series(C,G); 
Gpicl = feedback(Gpiol, 1), 

h2 = stepplot(Gpicl), 
p2 = getoptions(h2),
p2.Title.String = 'Unit Feedback Closed Loop Step Plot of PI-Controlled Soft Robot';
p2.TimeUnits= 'seconds',
setoptions(h2, p2),
stepinfo(Gpicl)

%% PD Secondary Loop feedforward
fprintf('Press any key to continue'),
pause,
clf reset,

Gpid =  pid(3.4993, 0.054765, 55.8988)   % pd controller with first order derivative

% use sisotool here to approximate inherent delay in TF model
Gnew = series(Gpid, Gpicl),
Gall = 1.171*feedback(Gnew,1);
%Pndall = pade(Gall, 2), 
h4 = stepplot(Gall),
p4 = getoptions(h4);
p4.Title.String = 'Cascade PID and PI Closed-Loop Step Response'; 
p4.TimeUnits = 'seconds';
setoptions(h4, p4); stepinfo(Gall)

fprintf('Press any key to continue'),
pause,
clf reset,
h7 = bodeplot(Gnew),
p7 = getoptions(h7),
p7.Title.String = 'Bode Plot of PID-PDControl Network',
setoptions(h7, p7)

%% Closed Loop Transfer Function
[num, den] = ss2tf(Pndall.a, Pndall.b, Pndall.c, Pndall.d)
[num2, den2] = ss2tf(Gall.a, Gall.b, Gall.c, Gall.d)
Model2 = tf(num2, den2)
%{
zpk =                0.017407 (s+2.222) (s-1.714) (s-1) (s+0.00866) (s^2 - 3s + 3)
  ----------------------------------------------------------------------------------------
  (s+3.232) (s-1.001) (s-1) (s+0.008617) (s^2 + 0.05945s + 0.02461) (s^2 + 3.055s + 3.141)

Model =       0.01741 s^6 - 0.06062 s^5 + 0.002221 s^4 + 0.2661 s^3 - 0.422 s^2 + 0.1952 s + 0.001722
          -----------------------------------------------------------------------------------------------
          s^8 + 4.355 s^7 + 1.752 s^6 - 9.396 s^5 - 7.903 s^4 + 9.429 s^3 + 0.5069 s^2 + 0.2538 s + 0.002156
%}
%% Discretize system
Ts = 0.05;
Del = exp(-2*s);
Pdel = pade(Del, 2),                   % approximate delay with 2nd order pade approximant 
Gdis = c2d(Ge, Ts, 'zoh'),              % Discrete Plant
deldis = c2d(Pdel, Ts, 'zoh')          % Discretized delay
Cdis = c2d(C, Ts, 'zoh'),              % Discrete PI controller
Cdis2 = pid(-5, -0.0433, 0, 0.2 + Ts/2, Ts)% Discrete PI controller, check
Gleaddis  = c2d(Glead, Ts, 'zoh'),     % Discretize lead controller
Gpddis =  pid(0.89238,0,9.7202,0.2 + Ts/2,Ts)          % Discrete pd controller
%% Root Locus
clear all; s = tf('s'); %s = -2.2222 + j*1.5159;
% PI Controller with Plant
Gpiol =   (-0.003*s^2 + 0.005115*s + 4.452e-05)* exp(-2*s)/...
                 (s^3 + 0.1128*s^2 + 0.001028*s);     % Open-loop TF
%s = tf('s');
Theta = (180/pi)*angle(Gpiol), M = abs(Gpiol); K = 1/M
% theta = -124.8829, pole = 3.2447
% new open-loop transfer function
Glpiol = (-0.003*s^2 + 0.005115*s + 4.452e-05)* exp(-2*s)/...
                 ((s^3 + 0.1128*s^2 + 0.001028*s)*(s + 3.2447)); 

%evaluate absolute value of new OLTF
Alpha = (180/pi)*angle(Glpiol), Mnew = abs(Glpiol); Knew = 1/M;

% we found K_lead = 6.5032 and angle(Glpiol) = 179.1175 deg
% new OLTF 
Gleadpi = (-0.0195*s^2 + 0.005115*s + 4.452e-05)* exp(-2*s)/...
                 ((s^3 + 0.1128*s^2 + 0.001028*s)*(s + 3.2447));


             
             
newSYS = 56.8557*((exp(-2*s)) *(s+2.2222)* (0.0006*s - 0.001028))/((s^2 + 0.1128*s + 0.001028)*(s+0.03948));
% Approximate delay in TF using the "pade" function
Pnd2 = pade(newSYS,2)    % second-order delay
% Plot root locus of new OL TF
h = rlocusplot(Pnd2);
p = getoptions(h);                                              % Get options for plot.
p.Title.String = 'Root Locus Plot For Lead Compensated System'; % Change title in options.
setoptions(h,p);          % Apply options to plot.
stepinfo(newSYS); 
% Generate step response of lead compensated system
subplot(211), step(Pnd2)
subplot(212), step(newSYS)
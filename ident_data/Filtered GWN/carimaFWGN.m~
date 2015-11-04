%% Identification of SRS Model (Summer/Fall 2015)
% © Olalekan Ogunmolu, 2015
 clear all; clc
cd('/home/lex/Documents/Matlab_Files/ident_data/Filtered GWN')

format short

input = load('/home/lex/Documents/Matlab_Files/ident_data/Filtered GWN/inputUGWN3.csv');
output = load('/home/lex/Documents/Matlab_Files/ident_data/Filtered GWN/outputUGWN3.csv');
powspec = load('/home/lex/Documents/Matlab_Files/ident_data/Filtered GWN/powspec.csv');
crest = load('/home/lex/Documents/Matlab_Files/ident_data/Filtered GWN/crest.csv');
autoc = load('/home/lex/Documents/Matlab_Files/ident_data/Filtered GWN/acf.csv');

Ts = 48/1000;

srs = iddata(output(1:length(input)), input, Ts);

present(srs);

%give names to the input and output channels and Time units

srs.InputName  = 'Current';
srs.OutputName = 'Fused Measurement';
srs.TimeUnit   = 'seconds';
srs.InputUnit  = 'mA';
srs.OutputUnit = 'mm';
srs.ExperimentName = 'Linear Filtered Gaussian White Noise Expt';

% split dataset in 60:40 ratio for estimation
[m,n] = size(srs)
srstest = srs(1:0.6*m);
srstrain = srs(1:0.4*m);
%% Input Signal Properties
format compact;
disp('Let''s see the power spectrum of the excitation signal for a minute');
plot(powspec),
legend('Power Spectral Density from LabVIEW''s FFT');

fprintf('\nPress any key to continue\n')
pause()

plot(crest),
legend('Crest of filtered white Gaussian noise signal')

fprintf('\nPress any key to continue\n')
pause(),

plot(autoc),
legend('autocorrelation of input excitation signal')


%% 
clf
time = srs.Ts * m;              %actual time for id expt
ptest = plot(srstest(4000:4200)), 
%xlim([0 203]); % 692.1 692.6]),
legend('Testing Data TIme Series');


fprintf('\nPress any key to plot Validation Time Series\n\n');
pause

clf
ptrain = plot(srstrain(600:1100)),
%xlim([0 53]) ;%axis([0 200 697.1 697.9]),
legend('Validation Time Series')

% data is not zero mean. So remove constant levels and make the data zero mean.
srs = detrend(srs);

srstest = srs(1:0.6*m);
srstrain = srs(1:0.4*m);

plot(srstest(4000:4500));
legend('Testing Data TIme Series');

fprintf('\nPress any key to plot Validation Time Series\n\n');
pause

plot(srstrain(1000:2000)),
legend('Validation Time Series')


%% Estimating Nonparametric Models
clf
%non-param fir model
srsi = impulseest(srs,[],'negative', impulseestOptions('RegulKernel','SE'));   
srsii = impulseest(srs);

% Show the 99.7% confidence bands
showConfidence(impulseplot(srsi), 3)

fprintf('\nPress spacebar key to continue\n \nEstimate impulse response of data\n'),
pause, clf;

showConfidence(impulseplot(srsii), 3)
grid on

% there is a 70-sample delay (dead-time) before the output responds to
% input or an 18-sample delay if we use the negative impulse function.

% there is probably a high-degree of feedback in the data. Future outputs
% are possibly coupled w/past inputs. We estimate the delay.
delay = delayest(srs)


%and find the probability of feedback with 
 
fprintf('\nPress spacebar key to continue\n'),
pause,

fprintf('\nThe feedback in the system is\n'),

fdbk = feedback(srs)
%a 100% feedback tells us there is high feedback in the data
% this corresponds to a highly unreliable kinect sensor used for the
% experiment


%The sample time of the data is 0.0556 second, while the process time
% constants might be much slower. We may detect some rather high
% frequencies in the output. In order to affirm this, let us first compute
% the input and output spectral analysis estimate: 
sy = spa(srs(:,1,[]));
su = spa(srs(:,[],1));

 
fprintf('\nPress spacebar key to continue\n \nCompute Spectrum \n'),
pause,


clf
wmin = 0.441; wmax = 4.53;   %choose desirable range for spectral plot
spectrum(sy,su)%, {wmin, wmax});
legend({'Output','Input'})
grid on
%the input has very little relative energy above 2.66 rad/sec while
% the output contains relatively larger values above that frequency. There
% seem to be some high frequency disturbances that may cause some problem for
% the model building 

%checkout input-output periodogram
clf
subplot(211)
periodogram(srs.y), 
legend('Output')

subplot(212)
periodogram(srs.u)
legend({'Input'})
grid on


fprintf('\nPress spacebar key to continue\n'),
pause,
fprintf('\n\nShow bodeplot of spectral analysis with 99.7 conf. bands\n')

clf
sz = spa(srs);
showConfidence(bodeplot(sz), 3);
grid on
%it seems HF disturbance is uncertain. We may do well to limit frequency
%range to 3.2 rad/s (= 10^.5 Hz)

%% Parametric Model for Process Behavior : SS Model
Nx = 1:10;
%I want the model order to be picked based on which Hankel matrix of the 
%ss has the appropriate log singular value
mdlss = ssest(srstest, Nx);  

%Fit to estimation data: 96.56% (prediction focus) FPE: 0.01198, MSE: 0.01196   
%--> Model is good
 
fprintf('\nPress spacebar key to continue\n'),
pause,

disp('\nDisplay state space model xtics\n');

mdlss,
 
fprintf('\nPress spacebar key to continue\n'),
pause,

disp('\nCompare state space model with estimation data\n');
pause(0.2),

compare(srstrain, mdlss); %, inf, 'Samples', 70:140);  % poor fit 25.38%


% check residue that model did not pick up
mode = 'CORR';
lags  = 30;

clf
resid(mdlss, srstrain, mode, lags)

fprintf('\nShow residue from input to residue using impulse response\n')
fprintf('\n Press any key \n');
pause

clf;
resid(mdlss, srstrain, 'IR');
legend('Residuals between Model and Data using impulse response model')

fprintf('\nShow residue from input to residue using fir response\n')
fprintf('\n Press any key \n');
pause

clf;
g = resid(mdlss, srstrain, 'FR');
legend('Residuals between Model and Data using freq response')
%% ARMAX Model
%I'll try with a 2nd order armax model for the test data and see if I can get
%good dynamics
opt = armaxOptions;
opt.Focus = 'stability';                  %estimate with prediction method and enforce model stability
opt.InitialCondition = 'estimate';        %estimate init condo with best least sq. fit
opt.SearchOption.MaxIter = 50;
opt.SearchOption.Tolerance = 1e-5;
opt.EstCovar = true;
opt.Display   = 'on';
opt.Regularization.Lambda = 100;
opt.Regularization.R = [1000];
opt.Regularization.Nominal = 'zero';

Nu = [10];                     %input offset to ignore, 
Ne = [1];                     %Number of experiments
Ny = [10];                     %Output offset to ignore

opt.InputOffset = [Nu];

na=2 ;  nb=2  ; nc=2  ; nk=1;                            %we know delay is 1 sec
mdlarmax = armax(srstest, [na,  nb ,  nc,   nk], opt);   %Fit to estimation data: 96.56%, FPE: 0.0119926

fprintf('\nPress spacebar key to continue\n'),
pause(1),
disp(mdlarmax);


disp('\nFit to Estimation Data: '), mdlarmax.Report.Fit,

pause(2),

clf
fprintf('\nPress spacebar key to continue\n'),
pause(3),

%check if model rates well with validation data
compare(srstrain, mdlarmax)                 %fit is 22.13%

pause(5)

resid(mdlarmax, srstrain, 'IR')
legend('impulse response from the input to the residuals')

pause(5)

resid(mdlarmax, srstrain, 'corr')
legend('Correlation analysis performed')

pause(5)

%frequency response from input to residue
resid(mdlarmax, srstrain, 'fr')
legend('frequency response from input to residuals')
%% Clearly the high-frequency noise in the output is an issue. Let's try
%decimating the data by a factor of 3

if exist('resample', 'file')==2
    %use resample comand from signal processing toolbox
    srsd = resample(srs, 4, 12);
else
    %use slower alternatiove
    srsd = idresample(srs, 4, 12);
end

[m,n] = size(srsd);
srsdtest = srsd(1:0.6*m);
srsdtrain = srsd(1:0.4*m);

% Trying to find a good structure for estimated data
Imp2 = impulseest(srsdtest);
showConfidence(impulseplot(Imp2, 30), 3);

fprintf('\nPress spacebar key to continue\n'),
pause(1),

delay2 = delayest(srsdtest)
%delay is now 27 samples

fprintf('\nPress spacebar key to continue\n'),
pause(1),

fdbl2 = feedback(srsdtest)                %still effin hundred

fprintf('\nPress spacebar key to continue\n'),
pause(1),

mdlarmax2 = armax(srsdtest, [na, nb, nc, nk])

pause(2)

mdlarmax2.Report.Fit

%compare fidelity of new model with previous models
M = inf;
compare(srstest, mdlss, mdlarmax, mdlarmax2, M)

%% Generate ARIMAX model for the slowly varying disturbance in system
mdlarimax = armax(srstest, [na,  nb ,  nc,   nk], 'IntegrateNoise', true, opt);
%Fit to estimation data: 96.56%, FPE: 0.0119994

compare(mdlarimax, srstrain, mdlarmax, M);

%check if model rates well with validation data
compare(srstrain, mdlarimax)                 %fit is 22.13%

resid(mdlarimax, srstrain, 'IR')
legend('impulse response from the input to the residuals')

pause(2)

resid(mdlarimax, srstrain, 'corr')
legend('Correlation analysis performed')

%% Find transients
clc, clf
subplot(211)
stepplot(mdlarmax),
legend('ARMAX Model', 'location', 'best')

subplot(212)
stepplot(mdlarimax)
legend('ARIMAX Model', 'location', 'best')

pause(2)

clf;
fprintf('\n\nShow Pole-Zero Map for mdlarmax and mdlarimax')

iopzmap(mdlarmax)
legend('Pole-Zero Map for mdlarmax')

pause(5)

clf
iopzplot(mdlarimax)
legend('Pole-Zero Map for mdlarimax')

pause(5)
clf;
fprintf('\n\nShow Pole-zero map for ss estimated model')
iopzmap(mdlss)
legend('Pole-Zero Map for State Space armax Model')

%% Check out O-L Bode Plots
clf, clc
p_ss = bodeplot(mdlss), 
legend('Bode plot for open-loop ss estimated model', 'location', 'best')

pause(4)

disp('Show 95% confidence intervals')

showConfidence(p_ss)

% Repeat for mdlarmax 
clf, clc
p_armax = bodeplot(mdlarmax), 
legend('Bode plot for open-loop armax estimated model', 'location', 'best')

pause(4)

disp('Show 95% confidence intervals')

showConfidence(p_armax, 50)

% Repeat for mdlarimax 
clf, clc
p_arimax = bodeplot(mdlarimax), 
legend('Bode plot for open-loop arimax estimated model', 'location', 'best')

pause(4)

disp('Show 95% confidence intervals')

showConfidence(p_arimax, 50)

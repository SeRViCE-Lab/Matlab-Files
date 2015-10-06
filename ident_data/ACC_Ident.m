%% Data Preprocessing
%
clc, clf;
format compact
cd('/home/lex/Documents/Matlab_Files/ident_data/Sawtooth Waveform')

input = load('currentval.csv');
output = load('fusion.csv');

save accconfdata.mat output input;

load accconfdata.mat;

dlmwrite('accdata.dat', [output, input])

load accdata.dat;
% 
% Display the raw data;
%
fprintf('\n');
fprintf('Step 1: Display the raw data (output/input)\n\n');
clf reset,
subplot(2,1,1),plot(accdata(:,1), 'LineWidth', 2.5)
subplot(2,1,1),title('Fused Measurement Signal'), hold on
xlabel('Samples'), ylabel('Measurement(mm)')
xlim([0 10128]), %ylabel([684, 696])
hold off

subplot(2,1,2),plot(accdata(:,2), 'LineWidth', 1.2), hold on
title('Input Current (mA)'), 
xlabel('Samples'), ylabel('Current(mA)'),
xlim([0 10128])
hold off
fprintf('        [Press ANY key to continue...]\n\n');
pause,
%
% Remove data means;
%
fprintf('Step 2: Remove data means and then data trend\n\n');
fprintf('   2.1  Remove data means\n\n');
gas_mean_removed=detrend(accdata,0);
clf reset,
subplot(2,1,1),plot(gas_mean_removed(:,1), 'LineWidth', 2.5), hold on
title('Mean Removed Output'), 
xlabel('Samples'), ylabel('Measurement(mm)')
xlim([0 10128]), %ylim([684, 696]), hold off

subplot(2,1,2),plot(gas_mean_removed(:,2), 'LineWidth', 1.5), hold on
title('Mean Removed Input'), 
xlabel('Samples'), ylabel('Current (mA)'),
xlim([0 10128])
hold off
fprintf('        [Press ANY key to continue...]\n\n');
pause,
%
% Remove data trends;
%
fprintf('   2.2  Detrended data\n\n');
gas_detrended=detrend(accdata);
clf reset,
subplot(2,1,1),plot(gas_detrended(:,1), 'LineWidth', 2.5), hold on
title('Detrended Output Signal'),
xlabel('Samples'), ylabel('Measurement(mm)')
xlim([0 10128])

subplot(2,1,2),plot(gas_detrended(:,2), 'LineWidth', 1.2), hold on
title('Detrended Input'), 
xlabel('Samples'), ylabel('Current (mA)'),
xlim([0 10128])
fprintf('        [Press ANY key to continue...]\n\n');
pause,
%
% Display the cross correlation function between 
% the raw input an doutput data;
%
fprintf('Step 3: Correlation analysis for raw input/output\n\n');
input_length=length(accdata);
z1=gas_detrended(:,1);
z2=gas_detrended(:,2);
vrawout=input_length*var(accdata(:,1));
vrawin=input_length*var(accdata(:,2));
varz1=input_length*var(z1);
varz2=input_length*var(z2);
ccf_raw=xcorr(accdata(:,1)-mean(accdata(:,1)),accdata(:,2)-mean(accdata(:,2)));
ccf_detrended=xcorr(z1-mean(z1),z2-mean(z2));
ccf_raw=ccf_raw/sqrt(vrawout*vrawin);
ccf_detrended=ccf_detrended/sqrt(varz1*varz2);

clf reset,
subplot(211),
plot(linspace(1, length(ccf_raw), length(ccf_raw) ), ccf_raw, 'LineWidth', 1.2),
 hold on, title('CCF [Raw Input-Output]');
xlabel('Lags'), grid,
xlim([0,25000]), hold off

subplot(212), plot(linspace(1, length(ccf_detrended), length(ccf_detrended) ),...
      ccf_detrended, 'LineWidth', 1.2), hold on
title('CCF [Detrended Input-Output]'),
xlabel('Lags'), grid, hold off
fprintf('        [Press ANY key to continue...]\n\n');
pause,

%
% Display the auto-correlation function the detrended input 
% and the prewhitended input;
%
fprintf('Step 4: Correlation analysis for detrnded\n');
fprintf('        and pre-whitened input and output\n\n');
% 
% Display the auto-correlation function of the detrended input; 
%{
fprintf('   4.1  Auto-corr analysis of detrended input\n\n');
acf_in=xcorr_acse(z2,z2);
subplot(111),
plot(acf_in),
xlabel('lags'),
title('Auto-Correlation Function of Detrended Input'),
fprintf('        [Press ANY key to continue...]\n\n');
pause,
%}
%xlabel('lag')
%
% Display the auto-correlation function of the pre-whitened input; 
%
reply=1;
while reply~=0
fprintf('   4.2  Auto-corr analysis of prewhitened input\n\n');
  fprintf('   (a)  Select model order for input prewhitening:\n'); 
  m=input('        m=');
  th=ar(z2,m);
  model_ar=ar(z2,m);
  fprintf('\n' );
  fprintf('   The filter model for input prewhitening:\n\n');
  fprintf('   '), model_ar,
  fprintf('\n');
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
  fprintf('   (b)  Prewhitened output and input\n\n');
  b=[0 th(3,:)];
  prewhitened_opt=z1+filter(b,1,z1);
  e=pe(model_ar,z2);
  clf reset,
  subplot(211),plot(prewhitened_opt),
  subplot(211),title('Pre-whitened Output #1'),
  subplot(212),plot(e),
  subplot(212),title('Pre-whitened Input #1'),
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
  fprintf('   (c)  Input whiteness test\n\n');
  subplot(111),
  [e,r]=resid_acse(z2,th);
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
%  fprintf('        Display the whiteness test again!\n\n');
%  subplot(111),
%  [e,r]=resid(model_ar,z2);
%  fprintf('        [Press ANY key to continue...]\n\n');
%  pause,
%
% Display the cross correlation function the pre-whitened
% input and output; 
%
  fprintf('   4.3  Cross correlation function between\n');
  fprintf('        the prewhitened input and output\n\n');
  ccf_opt_ipt=crosscorr_acse(prewhitened_opt,e,'coeff');
  x1=ccf_opt_ipt(296:320,1);
  time=[0:24];
  subplot(111),
  plot(time,x1,'-',time,x1,'o')
  title('CCF Between the Prewhitened Input and Output'),
  xlabel('lags')
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
  clear x1
  fprintf('        Are you satisfied with the results?\n');
  fprintf('        Please type a number to confirm!\n\n');
  fprintf('        0--YES; Other Numbers--NO!\n\n');
  reply=input('        Your choice is [0/others]: ');
  fprintf('\n\n');
  if reply~=0
     fprintf('        *****  Repeat Step 4.2! *****\n\n');
     reply=1;
  end
end
%
% Cross Spectral Analysis 
%  
fprintf('Step 5: Cross spectral analysis for\n');
fprintf('        detrended input and output\n\n');
cross_in=[z2;z2;z2];  % Periodically extension for input;
cross_out=[z1;z1;z1]; % Periodically extension for output;
subplot(111),
spectrum_acse(cross_in,cross_out);
fprintf('        [Press ANY key to continue...]\n\n');
pause,
fprintf(' \n');
%
% System Identification 
%  
fprintf('Step 6: System identification\n\n');
reply=1;
while reply~=0
  fprintf('   6.1  Model order selection\n\n');
   
  no=input('        Model order for output, no=');
  ni=input('        Model order for input, ni=');
  nn=input('        Model order for noise, nn=');
  nd=input('        Time delay for input, nd=');
  fprintf('\n\n');
  order=[no ni nn nd];
  th=armax_acse(gas_detrended,order);
  model_armax=armax(gas_detrended,order);
  fprintf('\n');
  fprintf('   The armax model for prewhitened input-output:\n\n')
  fprintf('   '),model_armax
  fprintf('\n');
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
  e=resid_acse(gas_detrended,th);
  fprintf('   6.2  Model reponses\n\n');
%
% Predicted output
%
  fprintf('   (a)  One-step-ahead predicted outputs\n\n');
  predicted_opt=z1-e;
  time=[1:1:length(accdata)]';
  subplot(111),
  plot(time,z1,'-',time,predicted_opt,'*')
  title('Measurement("-"); One-Step-Ahead Predicted Output("*")'),
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
%
% Deterministic output
%
  fprintf('   (b)  Model predicted outputs\n\n');
  deterministic_opt=idsim_acse(z2,th);
  subplot(111),
  plot(time,z1,'-',time,deterministic_opt,'+'),
  title('Measurement("-"); Model Predicted Output("+")'),
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
%
  fprintf('   (c)  Predicted residuals and errors\n');
  de=z1-deterministic_opt;
  time=[1:1:length(accdata)]';
  clf reset,
  subplot(211),plot(e)
  subplot(211),title('One-Step-Ahead Prediction Residuals'),
  subplot(212),plot(de)
  subplot(212),title('Model Prediction Errors'),
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
  fprintf('   (d)  Impulse response\n\n');
  z3=[4.0,zeros(1,length(accdata)-1)]';
  impulse_response=idsim_acse(z3,th);
  impy=[impulse_response(1:25,1)];
  time=[0:1:24]';
  subplot(111),
  plot(time,impy,'+',time,ccf_opt_ipt(296:320,1),'-', ...
                     time,ccf_opt_ipt(296:320,1),'o'),
  title('Impulse Response "+" and CCF "-" Between Prewhitened Output and Input'),
  fprintf('        [Press ANY key to continue...]\n\n');
  pause,
  fprintf('   6.3  Model validity test\n\n');
  subplot(111),
  [e,r]=resid_acse(gas_detrended,th);
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
  reply=input('        Your choice is [0/others]: ');
  fprintf('\n\n');
  if reply~=0
     fprintf('        ***** Repeat Step 6.1! *****\n\n');
     reply=1;
  end
end
close;
% HEAT_EXCHANGER_MODEL
% This programme provides an example to fit an ARX or ARMAX
% model for a heat exchanger model based on the band-limied
% white noise driven input-output data.
% 
% One-step-ahead prediction compared with the true measurements
% are plotted.
%
% The parameters for output regression terms, input regression
% terms and moving average terms are stored in A and B, respectively.
%
% This programme can be used for fitting general ARX and ARMAX 
% models by loading the associated data files.
%
% See also ARX, ARMAX, ARX_ACSE, ARMAX_ACSE, THETA for further help.
%
clc;
clear;
load blwn_driven_in;
load blwn_driven_out;
%
u=blwn_driven_in(2,:);
y=blwn_driven_out(2,:);
u1=u(1:500);
u2=u(501:1000);
y1=y(1:500);
y2=y(501:1000);
%
if size(u,1)==1
   u1=u1';
   u2=u2';
end
if size(y,1)==1
   y1=y1';
   y2=y2';
end
ny1=length(y1);
ny2=length(y2);
z1=detrend([y1 u1]);
z2=detrend([y2 u2]);
%
fprintf('\n\n');
repeat=1;
while repeat~=0
  fprintf('(1) Deterimine the model structure: ARX/ARMAX\n\n');
  structure=0;
  while (structure~=1) & (structure~=2) 
    fprintf('    1 -- ARX; 2 -- ARMAX\n\n');
    structure=input('    Your choice is [1/2]: ');
    fprintf('\n\n'); 
    if structure <1 | structure >2
      fprintf('    You should type 1 or 2, try again!\n\n');   
      structure=0;
    end
  end
%
  fprintf('(2) Determining the model orders:\n\n');
  if structure==1
     na=input('    Model order for output, na=');
     nb=input('    Model order for input, nb=');
     nk=input('    Time delay, nk=');
     order=[na nb nk];
     th=arx_acse(z1,order);
     model=arx(z1,order);
     fprintf('\n');
     fprintf('   The ARX model identified is as follows:\n\n');
     fprintf('   '),model
     fprintf('\n');
     fprintf('        [Press ANY to continue...]\n\n\n');
     pause,
  elseif structure==2
     na=input('    Model order for output, na=');
     nb=input('    Model order for input, nb=');
     nc=input('    Model order for noise, nc=');
     nk=input('    Time delay, nk=');
     order=[na nb nc nk];
     th=armax_acse(z1,order);
     model=armax(z1,order);
     fprintf('\n');
     fprintf('   The ARMAX model identified is as follows:\n\n');
     fprintf('   '),model
     fprintf('\n');
     fprintf('        [Press ANY to continue...]\n\n\n');
     pause,
  end
  [A,B,C,D]=polydata(model);
  fprintf('(3) Model validity test:\n\n');
  clf reset,
  [e,r]=resid_acse(z1,th);
  fprintf('        [Press ANY to continue...]\n\n\n');
  pause,
%
  fprintf('(4) How many steps do you want to predict?\n\n');
  np=input('    np-step-ahead prediction, np='); 
  fprintf('\n\n\');
%
  fprintf('(5) %u-step-ahead-prediction:\n\n',np);
  yp=predict(model,[y2,u2],np);
  err=y2-yp;
%
  clf reset,
  plot(ny1+1:ny1+ny2,y2(1:ny2),ny1+1:ny1+ny2,yp(1:ny2),'--');
  title('Multi-Step-Ahead Prediction (the Green Dashed Line)'), 
  fprintf('    [Press ANY key to continue...]\n\n');
  pause,
  fprintf('    %u-step-ahead-prediction error:\n\n',np);
  clf reset,
  plot(ny1+1:ny1+ny2,err(1:ny2)),
  title('Multi-Step-Ahead Prediction Errors'), 
  fprintf('    [Press ANY key to continue...]\n\n');
  pause,
  fprintf('    STD of %u-step-ahead-prediction errs: %f\n\n',...
               np,std(err));
  fprintf('\n');
  fprintf('(6) QUIT or REPEAT this programme ?\n\n');
  fprintf('    ***** 0--QUIT *****\n');
  fprintf('    ***** Other Numbers--REPEAT *****\n\n');
  repeat=input('    You choice is [0/others] : ');
  fprintf('\n\n');
  if repeat~=0
    fprintf('     Repeat the modelling process.\n\n\n\n');
  end    
end
close,


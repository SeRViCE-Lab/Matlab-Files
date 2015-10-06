% SUNSPOT_ARMA
% This programme provides an example to fit an AR model
% for the sunspot time series. This includes:
% 
% (1) Fit an ARMA model;
% (2) One-step-ahead prediction;
% (3) Multi-step-ahead prediction.
%
% The parameters for AR and MA regression terms are
% stored in A and C, respectively.
%
% This programme can be used for fitting general ARX models 
% by changing the name of the data file.
%
% Local functions:
%    RESID_ACSE
%
% See AR, ARX, ARMAX, AR_ACSE, ARX_ACSE, ARMAX_ACSE, THETA
% for further help.
%
clc;
clear;
load sunspot_yearly.dat,
load sunspot_monthly.dat,
%
data_index=0;
while (data_index~=1) & (data_index~=2)
   fprintf('\n\n');
   fprintf('(1) Please select the data file:\n\n');
   fprintf('    1-Yearly recorded sunspot time series.\n');
   fprintf('    2-Monthly recorded sunspot time series.\n\n');
   data_index=input('    Your selection is : ');
   if data_index==1
      y=sunspot_yearly(1:300,2);
   elseif data_index==2
      y=sunspot_monthly(1:3000,3);
   else
      fprintf('\n');
      fprintf('     You should type 1 or 2, try again!\n\n')
   end
end
if size(y,1)==1
   y=y';
end
ny=length(y);
fprintf('\n');
%
repeat=1;
while repeat~=0
   fprintf('(2) How many steps do you want to predict?\n\n');
   np=input('    np-step-ahead prediction, np='); 
   fprintf('\n\n');
   fprintf('(3) Determine the moder order\n\n');
   na=input('    The order for AR, na=');
   nc=input('    The order for MA, nc=');
   fprintf('        \n');
   fprintf('   ******************************\n');
   fprintf('   *                            *\n'); 
   if na <10
      fprintf('   *  Process order : na=%u      *\n',na);
   else 
      fprintf('   *  Process order : na=%u     *\n',na);
   end   
   if nc<10
      fprintf('   *  Noise order :   nc=%u      *\n',nc);
   else
      fprintf('   *  Noise order :   nc=%u     *\n',nc);
   end
   if np<10
      fprintf('   *  %u-step-ahead prediction   *\n\',np);
   else
      fprintf('   *  %u-step-ahead prediction  *\n\',np);
   end   
   fprintf('   *                            *\n'); 
   fprintf('   ******************************\n\n');
   fprintf('\n');
   fprintf('(4) Parameter estimation:\n\n');
   if data_index==2 & na+nc>10
      fprintf('    The progamme is computing the parameters.\n');
      fprintf('    Please wait a moment ...\n\n'); 
   end
   th=armax_acse(y,[na,nc]);
   model=armax(y,[na,nc]);
   [A,B,C,D]=polydata(model);
   fprintf('    The AR model identified is as follows:\n\n');
   fprintf(' '),model
   fprintf('    [Press ANY key to continue...]\n\n');
   pause,
   fprintf(' \n');
%      
%  Multi-step-ahead prediction
%
   fprintf('(5) %u-step-ahead-prediction:\n\n',np);
   if data_index==2 & na+nc>30
      fprintf('    Preparing for multi-step-ahead prediction.\n'); 
      fprintf('    Please wait a while ...\n\n'); 
   end
   yp=predict(model,y,np);
   err=y-yp;
%
   clf reset,
   plot(1:ny,y(1:ny),1:ny,yp(1:ny),'--');
   title('np-Step-Ahead Prediction (the Green Dashed Line)'), 
   fprintf('    [Press ANY key to continue...]\n\n');
   pause,
   fprintf('    %u-step-ahead-prediction error:\n\n',np);
   clf reset,
   plot(1:ny,err(1:ny)),
   title('np-Step-Ahead Prediction Error'), 
   fprintf('    [Press ANY key to continue...]\n\n');
   pause,
   fprintf('    STD of %u-step-ahead-prediction err: %f\n\n',...
                np,std(err));
   fprintf('\n');
   fprintf('(6) Model validity test:\n\n');
   [e,r]=resid_acse(y,th);
   fprintf('    [Press ANY key to continue...]\n\n');
   pause,
%   fprintf('    Model validity test again:\n\n');
%   [e,r]=resid(model,y);
%   fprintf('    [Press ANY key to continue...]\n\n\n');
%   pause,
   fprintf('(7) QUIT or REPEAT this programme ?\n\n');
   fprintf('    ***** 0--QUIT *****\n');
   fprintf('    ***** Other Numbers--REPEAT *****\n\n');
   repeat=input('    You choice is [0/others] : ');
   fprintf('\n\n\n');
   if repeat~=0
      fprintf('     Repeat the modelling process.\n\n');
   end    
close,
end  
   


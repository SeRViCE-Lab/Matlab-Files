function [e,r]=resid_acse(z,th,M,maxsize)
%RESID_ACSE  
%   Computes and tests the residuals associated with a model.
%   E = RESID_ACSE(Z,TH)
%
%   Z : The output-input data Z=[y u], with y and u being column vectors.
%   For multi-variable systems Z=[y1 y2 .. yp u1 u2 ... un]. 
%   For time-series Z=y only.
%   TH: The model to be evaluated on the given data set. (For the theta 
%   format, see also THETA.)
%
%   E : The residuals associated with TH and Z. [RESID(Z,TH); just performs
%   and displays the tests, without returning any data.]
%
%   The autocorrelation function of E and the cross correlation between
%   E and the input(s) is computed and displayed. 99 % confidence limits 
%   for these values are also given (based on the hypothesis
%   that the residuals are white and independent of the inputs). These 
%   functions are given up to lag 25, which can be changed to M by
%   E = RESID(Z,TH,M). The correlation information can be saved and re-
%   plotted by
%   [E,R] = RESID_ACSE(Z,TH).         RESID_ACSE(R);
%
%   E = RESID_ACSE(Z,TH,M,MAXSIZE) changes the memory variable MAXSIZE from
%   its default value. See also AUXVAR.
%   See also PE.

%   Copyright (c) 1986-98 by The MathWorks, Inc.
%   $Revision: 2.4 $  $Date: 1997/12/02 03:42:42 $

if nargin < 1
   disp('Usage: RESID(DATA,MODEL)')
   disp('       E = RESID(DATA,MODEL,No_OF_LAGS,MAXSIZE)')
   return
end


% ** Set up the default values **
newplot;
[N,nz]=size(z);
maxsdef=idmsize(N);

if nargin<4, maxsize=maxsdef;end
if nargin<3, M=25;end
if maxsize<0,maxsize=maxsdef;end,if M<0, M=25;end
if nargin==1,[nz2,M1]=size(z);nz=sqrt(nz2);M=M1-2;N=z(1,M1-1);
ny=z(1,M1);nu=nz-ny;r=z(:,1:M);end
if nargin>1,
[N,nz]=size(z); nu=th(1,3);ny=nz-nu;
%
% *** Compute the residuals and the covariance functions ***
%
e=pe_acse(z,th);

if nz>1, e1=[e z(:,ny+1:nz)];else e1=e;end
r=covf(e1,M,maxsize);
end
nr=0:M-1;plotind=0;
for ky=1:ny
ind=ky+(ky-1)*nz;
%
% ** Confidence interval for the autocovariance function **
%
sdre=2.58*(r(ind,1))/sqrt(N)*ones(M,1);
if nz==1,spin=111;else spin=210+plotind+1;end,subplot(spin)
plot(nr,r(ind,:)'/r(ind,1),'-',nr,r(ind,:)'/r(ind,1),'o',nr,[ sdre -sdre]/r(ind,1),':r')
title(['Correlation Function of Residuals. Output # ',int2str(ky)])
%xlabel('lags')
plotind=rem(plotind+1,2);
if plotind==0,pause,newplot;end
end
nr=-M+1:M-1;
% *** Compute confidence lines for the cross-covariance functions **
for ky=1:ny
for ku=1:nu
ind1=ky+(ny+ku-1)*nz;ind2=ny+ku+(ky-1)*nz;indy=ky+(ky-1)*nz;
indu=(ny+ku)+(ny+ku-1)*nz;
   %corr 891007
   sdreu=2.58*sqrt(r(indy,1)*r(indu,1)+2*(r(indy,2:M)*r(indu,2:M)'))/sqrt(N)*ones(2*M-1,1); % corr 890927
   spin=210+plotind+1;subplot(spin)
   plot(nr,[r(ind1,M:-1:1) r(ind2,2:M) ]'/(sqrt(r(indy,1)*r(indu,1))),'-',...
   nr,[r(ind1,M:-1:1) r(ind2,2:M) ]'/(sqrt(r(indy,1)*r(indu,1))),'o', ...
   nr,[ sdreu -sdreu]/(sqrt(r(indy,1)*r(indu,1))),':r' )

   title(['Cross Corr. Function Between Input #' int2str(ku)...
         ' and Residuals From Output #',int2str(ky)])
   xlabel('lags')
   plotind=rem(plotind+1,2);if ky+ku<nz & plotind==0,pause,newplot;end
end,end
r(1,M+1:M+2)=[N,ny];
set(gcf,'NextPlot','replace');


function [th,ref]=ar_acse(y,n,approach,win,maxsize,T)
%AR_ACSE   Computes AR-models of signals using various approaches.
%   TH = AR_ACSE(y,n)  or  TH = AR_ACSE(y,n,approach)  or 
%   TH = AR_ACSE(y,n,approach,win)
%
%   TH: returned as the estimated parameters of the AR-model,see also THETA.
%
%   y: The time series to be modelled, a column vector
%   n: The order of the AR-model
%   approach: The method used, one of the following ones:
%      'fb' : The forward-backward approach (default)
%      'ls' : The Least Squares method
%      'yw' : The Yule-Walker method
%      'burg': Burg's method
%      'gl' : A geometric lattice method
%      For the two latter ones, reflection coefficients and loss functions
%      are returned in REFL by [TH,REFL]=AR(y,n,approach)
%      If any of these arguments end with a zero (like 'burg0'), the
%      computation of the covariance is suppressed.
%   win : Windows employed, one of the following ones:
%      'now' : No windowing (default, except when approach='yw')
%      'prw' : Prewindowing
%      'pow' : Postwindowing
%      'ppw' : pre- and post-windowing
%   for TH = AR_ACSE(y,n,approach,win,maxsize,T) see also AUXVAR.
%   See also IVAR and for the multi-output case ARX and N4SID.
%
%   See also ARX, ARMAX, ARX_ACSE, ARMAX_ACSE, THETA_ACSE

%   Copyright (c) 1986-98 by The MathWorks, Inc.
%   $Revision: 2.4 $  $Date: 1997/12/02 03:41:21 $

if nargin <2
   disp('Usage: TH = AR(Y,ORDER)')
   disp('       TH = AR(Y,ORDER,APPROACH,WINDOW,MAXSIZE,T)')
   disp('       APPROCH is one of ''fb'', ''ls'', ''yw'', ''burg'', ''gl''.')
   disp('       WINDOW is one of ''now'', ''prw'', ''pow'', ''ppw''.')
   return
end

% Some initial tests on the input arguments
[Ncap,ny]=size(y);Ncap0=Ncap;
maxsdef=idmsize(Ncap,n);
if nargin<6,T=[];end
if nargin<5,maxsize=[];end
if nargin<4,win=[];end
if nargin<3,approach=[];end
if isempty(T),T=1;end,if T<0,T=1;end
if isempty(maxsize), maxsize=maxsdef;end,if maxsize<0,maxsize=maxsdef;end
if isempty(win);win='now';end ,if win<0,win='now';end
if isempty(approach),approach='fb';end,if approach<0,approach='fb';end
pt=1; %pt=1 means that covariances should be computed
if approach(length(approach))=='0',pt=0;end
if length(approach)>1,appr=approach(1:2);else appr=[approach ' '];end
ap=[appr;appr;appr;appr;appr;appr;appr;appr;appr;appr];
ss=(ap==['fb';'FB';'yw';'YW';'ls';'LS';'bu';'BU';'gl';'GL']);
if sum(ss(:,1).*ss(:,2))==0,error('The input argument APPROACH must be one of the following (within quotes): ''fb'', ''ls'',''yw'', ''burg'', ''gl'''),end
ap=[win;win;win;win]; ss=(ap==['prw';'pow';'ppw';'now']);
if sum((ss(:,1).*ss(:,2)).*ss(:,3))==0,error('The input argument WIN must be one of the following (within quotes): ''now'',     ''prw'', ''pow'', ''ppw'''),end
[Ncap,ny]=size(y); if ny>Ncap, y=y.';end
if min(Ncap,ny)>1, error('Only single time series can be handled'),end

if appr=='yw' | appr=='YW',win='ppw';end
if win=='prw' | win=='ppw',y=[zeros(n,1);y];nstart=n+1;nend=length(y);end
if win=='pow' | win=='ppw',y=[y;zeros(n,1)];nend=length(y)-n;end

[Ncap,ny]=size(y);
nstart=1;nend=Ncap;
%build up the theta-matrix
%
I=[0 T 0 n 0 0];
if pt~=0,th=zeros(n+3,max(7,n));else th=zeros(3,max(7,n));end
th(1,1:6)=I;
ti=fix(clock);ti(1)=ti(1)/100;
th(2,2:6)=ti(1:5);th(2,7)=12;

% First the lattice based algorithms

if appr=='bu' | appr=='BU' | appr=='gl' | appr=='GL'
   ef=y;eb=y;
   rho(1)=y'*y/Ncap;
   for p=1:n
   nef=ef(p+1:Ncap)'*ef(p+1:Ncap); neb=eb(p:Ncap-1)'*eb(p:Ncap-1);
   if appr=='gl' | appr=='GL',den=sqrt(nef*neb);end
   if appr=='bu' | appr=='BU',den=(nef+neb)/2;end
   r(p)=(-eb(p:Ncap-1)'*ef(p+1:Ncap))/den;
   efold=ef;
   ef(2:Ncap)=ef(2:Ncap)+r(p)*eb(1:Ncap-1);
   A(p)=r(p);
   eb(2:Ncap)=eb(1:Ncap-1)+conj(r(p))*efold(2:Ncap);
   A(1:p-1)=A(1:p-1)+r(p)*conj(A(p-1:-1:1));
   rho(p+1)=rho(p)*(1-r(p)*r(p));
   end
th(3,1:n)=A;
e=pe(y(nstart:nend),th);lam=e'*e/(nend-nstart+1-n);
th(1,1)=lam; th(2,1)=lam*(1+n/Ncap0)/(1-n/Ncap0);
ref=[0 r;rho];
if pt==0, return,end
end
% Now compute the regression matrix

nmax=n;
M=floor(maxsize/n);
R=zeros(n);F=zeros(n,1);
fb=((appr=='fb') | (appr=='FB'));
if fb,RB=zeros(n);FBc=zeros(n,1);yb=conj(y(Ncap:-1:1));end
for k=nmax:M:Ncap-1
   jj=(k+1:min(Ncap,k+M));
   phi=zeros(length(jj),n);if fb,phib=zeros(length(jj),n);end
   for k=1:n,phi(:,k)=-y(jj-k);end
   if fb,for k=1:n,phib(:,k)=-yb(jj-k);end,end
   R=R+phi'*phi;F=F+phi'*y(jj);
   if fb,RB=RB+phib'*phib;FBc=FBc+phib'*yb(jj);end
end
P=inv(R);
if appr~='bu' & appr~='gl'
if ~fb
   if Ncap>M , t=R\F;else t=phi\y(jj);end
end
if fb
   if Ncap>M/2 , t=(R+RB)\(F+FBc);else t=[phi;phib]\[y(jj);yb(jj)];end
end
th(3,1:n)=t.';
e=pe_acse(y(nstart:nend),th);lam=e'*e/(nend-nstart+1-n);
th(1,1)=lam; th(2,1)=lam*(1+n/Ncap0)/(1-n/Ncap0);
end
if pt~=0,P=lam*P;th(4:n+3,1:n)=P;end
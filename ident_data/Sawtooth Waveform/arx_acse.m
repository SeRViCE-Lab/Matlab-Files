function th=arx_acse(z,nn,maxsize,Tsamp,p);
%ARX_ACSE    Computes LS-estimates of ARX-models.
%   TH = ARX_ACSE(Z,NN)
%
%   TH: returned as the estimated parameters of the ARX model
%   A(q) y(t) = B(q) u(t-nk) + e(t)
%   along with estimated covariances and structure information.
%   For the exact structure of TH see also THETA.
%
%   Z : the output-input data Z=[y u], with y and u as column vectors.
%   For multivariable systems Z=[y1 y2 .. yp u1 u2 .. um]. For time series
%   Z=y only.
%
%   NN: NN = [na nb nk], the orders and delays of the above model.
%   For multi-output systems, NN has as many rows as there are outputs
%   na is then an ny|ny matrix whose i-j entry gives the order of the
%   polynomial (in the delay operator) relating the j:th output to the
%   i:th output. Similarly nb and nk are ny|nu matrices. (ny:# of outputs,
%   nu:# of inputs). For a time series, NN=na only.
%   Some parameters associated with the algorithm are accessed by
%   TH = ARX_ACSE(Z,NN,maxsize,T)
%   See also AUXVAR for an explanation of these and their default values.
%   See also ARXSTRUC, SELSTRUC for how to estimate the order parameters.
%
%   See also AR, ARMAX, ARX_ACSE, ARMAX_ACSE, THETA_ACSE, BJ, IV4, N4SID, 
%            OE, PEM.
%

%   Copyright (c) 1986-98 by The MathWorks, Inc.
%   $Revision: 2.4 $  $Date: 1997/12/02 03:39:50 $

if nargin <2
   disp('Usage: TH = ARX(Z,ORDERS)')
   disp('       TH = ARX(Z,ORDERS,MAXSIZE,T)')
   return
end

% *** Set up default values **
[Ncap,nz]=size(z); [nr,nc]=size(nn);
maxsdef=idmsize(Ncap); % Mod 891007
if nargin<5, p=1;end
if nargin<4, Tsamp=[];end
if nargin<3, maxsize=[];end
if isempty(Tsamp),Tsamp=1;end,if Tsamp<0,Tsamp=1;end
if isempty(maxsize),maxsize=maxsdef;end,if maxsize<0,maxsize=maxsdef;end
if p<0,p=1;end
if any(any(nn<0))
   error('All orders must be non-negative.')
end

if nz>Ncap,error('Data should be organized in column vectors!'),return,end
if nr>1, eval('th=mvarx(z,nn,maxsize,Tsamp);'),return,end
nu=nz-1;
if length(nn)~=1+2*(nz-1)
   disp(' Incorrect number of orders specified!')
   disp('For an AR-model nn=na'),disp('For an ARX-model, nn=[na nb nk]')
   error(' ')
end

if nz>1, na=nn(1);nb=nn(2:1+nu);nk=nn(2+nu:1+2*nu);
else na=nn(1); nb=0;,nk=0;end
n=na+sum(nb);
if n==0,
   th=[];
   if p>0, disp('All orders are zero. No model returned.');end
   return
end

% *** construct regression matrix ***
 
nmax=max([na+1 nb+nk])-1;
M=floor(maxsize/n);
R=zeros(n);F=zeros(n,1);
for k=nmax:M:Ncap-1
      jj=(k+1:min(Ncap,k+M));
      phi=zeros(length(jj),n);
      for kl=1:na, phi(:,kl)=-z(jj-kl,1);end
      ss=na;
      for ku=1:nu
           for kl=1:nb(ku), phi(:,ss+kl)=z(jj-kl-nk(ku)+1,ku+1);end
           ss=ss+nb(ku);
      end
      if Ncap>M | p~=0, R=R+phi'*phi; F=F+phi'*z(jj,1);end
end
%
% *** compute estimate ***
%
if Ncap>M, th=R\F; else th=phi\z(jj,1);end
if p==0, return,end
%
% proceed to compute loss function and covariance matrix
%
t=th;, clear th;

%
% build up the theta-matrix
%
if nu>0
     I=[1 Tsamp nu na nb 0 0 zeros(1,nu) nk];
else I=[1 Tsamp 0 na 0 0];end

n=na+sum(nb);
th=zeros(n+3,max([length(I) n 7]));
th(1,1:length(I))=I;
ti=fix(clock); ti(1)=ti(1)/100;
th(2,2:6)=ti(1:5);
th(2,7)=1;
th(3,1:n)=t.';
e=pe_acse(z,th);V=e'*e/(Ncap-nmax);
th(4:n+3,1:n)=V*inv(R); % mod 911209
th(1,1)=V;
th(2,1)=V*(1+n/Ncap)/(1-n/Ncap); %Akaike's FPE


function [y,ysd]=idsim_acse(ue,th,inhib)
%IDSIM  Simulates a given dynamic system.
%   Y = IDSIM_ACSE(UE,TH)
%
%   TH: contains the parameters of the model in the format described by 
%   (see also) THETA.
%
%   UE: the input-noise data Z=[u e]. For multi-variable systems
%   u=[u1 u2 ... un e1 e2 .. ep], with ui and ei being column vectors. 
%   The number of noise sources should equal the number of outputs. If
%   the e-vector(matrix) is omitted, a noise-free simulation is obtained.
%   The noise contribution is scaled by the variance information con-
%   tained in TH.
%
%   With  [Y,YSD] = IDSIM_ACSE(UE,TH)  the estimated standard deviation of the
%   simulated output, YSD, is also computed. This is currently only 
%   available for non-statespace models.
%   See also IDINPUT, MS2TH and POLY2TH for simulation and model creation.
%   See also COMPARE and PREDICT for model evaluation.

%   Copyright (c) 1986-98 by The MathWorks, Inc.
%   $Revision: 2.5 $  $Date: 1997/12/03 23:15:54 $

y=[];ysd=[];
if nargin<3,inhib=0;end
if nargin<2
   disp('Usage: Y = IDSIM(UE,TH)')
   disp('       [Y,YSD] = IDSIM(UE,TH)')
   return
end

if isthss(th),
   if nargout>1&~inhib,
     disp(str2mat('No standard deviations will be computed for',...
                  'multi-output and state-space models.',...
                  'You could use THSS2TH for a conversion to MISO-models.')) 
   end
   eval('y=idsimss(ue,th);')
   return
end
[N,nue]=size(ue);Ncap=N;
if nue>N, error(' The data should be organized in columns.'),return,end

T=gett_acse(th);
if T<0
  if inhib
     th=thc2thd(th,-T);
  else
     disp(str2mat('This is a continuous time model of input-output type.',...;
         'Sample it (using thc2thd) before applying idsim.'));
     error(' ')
  end
end
[a,b,c,d,f]=th2poly_acse(th);
nu=th(1,3);
LAM=th(1,1);
y=zeros(N,1);
if nue < nu
   error(['The number of specified input signals is less than what the',...
         ' model requires.'])
end
for k=1:nu
   y=y + filter(b(k,:),f(k,:),ue(:,k));
end

if nue>nu, y=y+sqrt(LAM)*filter(c,d,ue(:,nu+1));end

y=filter([1],a,y);

if nargout>1 % now we compute the standard deviation of y:
   na=th(1,4); nb=th(1,5:4+nu); nc=th(1,5+nu);
   nd=th(1,6+nu); nf=th(1,7+nu:6+2*nu); nk=th(1,7+2*nu:6+3*nu);
   n=na+sum(nb)+sum(nf);
   nmax=max([na nb+nk-ones(1,nu)  nf]);

   % *** Prepare for gradients calculations ***

   if na>0, yf=filter(-1,a,y);end
   for k=1:nu
      gg=conv(a,f(k,:));
      uf(:,k)=filter(1,gg,ue(:,k));
      if nf(k)>0, wf(:,k)=filter(-b(k,:),f(k,:),uf(:,k));end
   end

% *** Compute the gradient PSI. ***


      psi=zeros(Ncap,n);%jj=nmax+1:Ncap;
      for kl=1:na, psi(kl+1:Ncap,kl)=yf(1:Ncap-kl);,end
      ss=na;ss1=na+sum(nb);
      for ku=1:nu
             for kl=1:nb(ku), psi(kl+nk(ku):Ncap,ss+kl)=...
                                 uf(1:Ncap-kl-nk(ku)+1,ku);end
             for kl=1:nf(ku), psi(kl+1:Ncap,ss1+kl)=wf(1:Ncap-kl,ku);end
             ss=ss+nb(ku);ss1=ss1+nf(ku);
      end

%*** THe covariance of y:

  [par,P]=th2par(th);
  actpar=[1:na+sum(nb),na+sum(nb)+nc+nd+1:length(par)];
  [nrp,ncp]=size(P);
  if nrp>=length(actpar);
      P=P(actpar,actpar);
      ysd=sqrt(sum(((psi*P).*psi)'))';
  else
      ysd=[];
  end
end

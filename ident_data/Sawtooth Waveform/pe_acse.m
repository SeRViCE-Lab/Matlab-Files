function [e,v,w,a,b,c,d,f]=pe_acse(z,th)
%PE_ACSE   Computes prediction errors.
%   E = PE_ACSE(Z,TH)
%
%   E : The prediction errors
%   Z : The output-input data Z=[y u]
%   TH: The model. Format as in HELP THETA
%
%   A more complete set of information associated with the model TH and the
%   data Z is obtained by
%
%   [E,V,W,A,B,C,D,F] = PE_ACSE(Z,TH)
%
%   Here A,B,C,D,and F are the polynomials associated with TH and
%   W = B/F u and V = A[y - W].

%   Copyright (c) 1986-98 by The MathWorks, Inc.
%   $Revision: 2.5 $  $Date: 1997/12/02 03:40:16 $

if nargin < 2
   disp('Usage: E = PE(Z,TH);')
   return
end
[Ncap,nz]=size(z);
nu=th(1,3);
if isthss(th),
   [nr,nc]=size(th);
   nd=th(1,5);
   sspmod=getmfth(th);

   [etapar,P,lambda]=th2par(th);
   arg=getargth(th);
   T=th(1,2);
   if any(th(2,8)==[2 3]),Tmod=-1;else Tmod=abs(T);end
   if length(etapar)==0,etapar=0;end
   [a,b,c,d,k,x0]=feval(sspmod,etapar,Tmod,arg);
   [ny,nx]=size(c);
   if nu+ny~=nz, 
      disp('Incorrect number of data columns specified.') 
      disp(['Should be equal to the sum of the number of inputs'])
      disp('and the number of outputs.')
      error(' ')
   end 
   x=ltitr(a-k*c,[k b-k*d],z,x0);
   e=z(:,1:th(1,4))-(c*x')';
   if ~isempty(d),e=e-(d*z(:,th(1,4)+1:th(1,4)+th(1,3))')';end
else

   if nu~=nz-1
      disp('Incorrect number of data columns specified.') 
      disp(['Should be equal to the sum of the number of inputs'])
      disp('and the number of outputs.')
      error(' ')
   end
   [a,b,c,d,f]=th2poly_acse(th);[nu,nb]=size(b);[nu,nf]=size(f);
   nc=length(c);nd=length(d);
   ni=max([length(a)+nd-2 nb+nd-2 nf+nc-2 1]);
   v=filter(a,1,z(:,1));
   for k=1:nu, 
      if ~isempty(b)
          w(:,k)=filter(b(k,:),f(k,:),z(:,k+1));v=v-w(:,k);
      end
   end
   e=pefilt_acse(d,c,v,zeros(1,ni));
end

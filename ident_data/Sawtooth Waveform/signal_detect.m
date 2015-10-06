%
% SIGNAL_DETECT
%
% This programme provides an example on the detection of
% a signal buried in noise using cross-correlation.
% 
% The deterministic signal is
%     u(t)=cos(100*pi*t)
% The smaple data will be obtained with the sampling frequency
% Fs=2000Hz and 1000 points will be recorded and used.
%
% Two types of noise are considered:
% (1) Normal distribution with zero mean and the std=2.5;
% (2) Uniform distribution on [-4,4].
%
% Local functions:  XCORR_ACSE
%
clc,
clear;
repeat=1;
while repeat~=0
   noise_type=0;
   while (noise_type~=1) & (noise_type~=2)
      fprintf('\n\n');
      fprintf('(1) Please select the noise_type:\n\n');
      fprintf('    1-Normal distribution N(0,xigmma^2)\n');
      fprintf('    2-Uniform distribution on [-a,a]\n');
      fprintf('    [xigmma=2.5 and a=4 are recommanded]\n\n');
      noise_type=input('    Your selection is [1 or 2] : ');
      if noise_type==1
      fprintf('\n\n');
         fprintf('(2) Input a value for parameter "xigmma":\n\n');
         xigmma=input('    XIGMMA=[2.5/others] : ');
         n=1;
         for t=0:1/2000:0.5
            u(n)=cos(100*pi*t);
            rn(n)=xigmma*randn;
            y(n)=u(n)+rn(n);
            z(n)=0;
            n=n+1;
         end
      elseif noise_type==2
      fprintf('\n\n');
         fprintf('(2) Input a value for parameter "a" :\n\n');
         a=input('    A=[4/others] : ');
         n=1;
         for t=0:1/2000:0.5
            u(n)=cos(100*pi*t);
            ru(n)=-a+2*a*rand;
            y(n)=u(n)+ru(n);
            z(n)=0;
            n=n+1;
         end
      else
         fprintf('\n');
         fprintf('     You should type 1 or 2. Try again!\n\n')
      end
   end
   fprintf('\n');
%
   N=min([length(0:1/2000:0.5),1000,length(u)]);
   Ruy=xcorr_acse(u,y);
%
   ymax=max(abs(y));
   rmax=max(abs(Ruy));
   clf reset,
   subplot(3,1,1),plot(0:N-1,u(1:N),0:N-1,z(1:N));
   title('Detecting A Signal Buried in Noise Using Cross Correlation'), 
   axis([0 N -1.5 1.5]);
   ylabel('u'),
   text(750,-1.25,'u=cos(100*pi*t)');
   text(960,1.25,'(a)');
%
   subplot(3,1,2),plot(0:N-1,y(1:N),0:N-1,z(1:N));
   axis([0 N -1.2*ymax 1.2*ymax]);
   text(750,-1.0*ymax,'y=u+noise');
   text(960,1.0*ymax,'(b)');
   ylabel('y'),
%
   subplot(3,1,3),
   plot(0:N-1,Ruy(length(u):length(u)+N-1),0:N-1,z(1:N));
   axis([0 N -1.5*rmax 1.5*rmax]);
   ylabel('CCFuy'),
   text(750,-1.25*rmax,'Cross Correlation');
   text(960,1.25*rmax,'(c)');
%
   fprintf('    [Press ANY key to continue...]\n\n\n');
   pause,
   fprintf('(3) Satisfied with the results?\n\n');
   fprintf('    0--YES;  Orther Numbers--NO\n\n');
   repeat=input('    You choice is [0/others]: ');
   fprintf('\n\n');
   if repeat~=0
       fprintf('    *****  Repeat  *****\n');
   end
end
close,
% FREQ_ANALYSIS
% This programme is used for
% (1) Calculating the frequency response for mdels
%     (3.1),(3.2b) and (3.3).
% (2) Plotting the Bode diagrams for models
%     (3.1), (3.2b)
%
clear;
clc;
m=1;
for f=0.0:0.001:0.6
   rp0(m)=freqresp([0.8],[900,420,43,1],f);          % Model (3.1)
   rp2(m)=freqresp([0.435 -2.332 2.924], ...
                   [3957.2,80,11.2],f);              % Model (3.2b)
   rp3(m)=freqresp([0.00024,-0.00059,-0.00261,0.00566], ...
                   [6.4323,2.9814,0.3044,0.0072],f); % Model (3.3)
m=m+1;
end
ff=0.0:0.001:0.6;
mf=length(ff);
fprintf('\n\n');
fprintf(' Freqency Response of the two models:\n');
fprintf(' (1) G(s)=0.8/[(30s+1)(10s+1)(3s+1)]\n');
fprintf(' (2) G(z)=(0.000138z+0.000575)/(z^2-1.950z+0.9785)\n\n');
%
clf reset,
plot(ff,rp0(1:mf),ff,rp2(1:mf),'--'),
title('Freq. Resp. of G(s)=0.8/[(30s+1)(10s+1)(3s+1)](the solid blue line)'),
xlabel('Frequency [Hz]'); grid,
fprintf('     [Press any key to contunue...]\n\n');
pause,
%
fprintf(' Freqency Response for the two models:\n');
fprintf(' (1) G(s)=0.8/[(30s+1)(10s+1)(3s+1)]\n');
fprintf(' (2) G(z)=0.001(0.1039z^2+0.4734z^1+0.1031)/(z^3-2.5886z^2+2.2166z-0.6271)\n\n');
clf reset,
plot(ff,rp0(1:mf),ff,rp3(1:mf),'--'),
title('Freq Resp of G(z^-^1)=10^-^3(0.104z^-^1+0.47z^-^2+0.1z^-^3)/(1-2.5886z^-^1+2.2166z^-^2-0.6271z^-^3)'),
xlabel('Frequency [Hz]'); grid,
fprintf('     [Press any key to contunue...]\n\n');
pause,
%
fprintf(' Bode Diagram for the model\n');
fprintf(' G(s)=0.8/[(30s+1)(10s+1)(3s+1)]\n\n');
clf reset,
bode([0.8],[900,420,43,1]),grid,
title('Bode Diagram of G(s)=0.8/[(30s+1)(10s+1)(3s+1)]'),
fprintf('     [Press any key to contunue...]\n\n');
pause,
%
fprintf(' Bode Diagram for the model\n');
fprintf(' G(z)=(0.000138z+0.000575)/(z^2-1.975z+0.9788)\n\n');
clf reset,
bode([0.435 -2.332 2.924],[3957.2,80,11.2]),grid,
title('Bode Diagram of G(z^-^1)=(0.000138+0.000575z^-^1)/(1-1.975z^-^1+0.9788z^-^2)'),
fprintf('     [Press any key to contunue...]\n\n');
pause,
close,

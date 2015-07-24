%% Filter examples
load bostemp
days = (1:31*24)/24;
plot(days, tempC), axis tight;
ylabel('Temp (\circC)');
xlabel('Time elapsed from Jan 1, 2011 (days)');
title('Logan Airport Dry Bulb Temperature (source: NOAA)');
%A Moving Average Filter
clf;
hoursPerDay = 24;
coeff24hMA = ones(1, hoursPerDay)/hoursPerDay;

avg24hTempC = filter(coeff24hMA, 1, tempC);
plot(days, [tempC avg24hTempC]);
legend('Hourly Temp','24 Hour Average (delayed)','location','best');
ylabel('Temp (\circC)');
xlabel('Time elapsed from Jan 1, 2011 (days)');
title('Logan Airport Dry Bulb Temperature (source: NOAA)');
% Model the delay
clf
fDelay = (length(coeff24hMA)-1)/2;
plot(days, tempC, 'b', ...
     days-fDelay/24, avg24hTempC, 'g');
axis tight;
legend('Hourly Temp','24 Hour Average','location','best');
ylabel('Temp (\circC)');
xlabel('Time elapsed from Jan 1, 2011 (days)');
title('Logan Airport Dry Bulb Temperature (source: NOAA)');
%Extract Avg Difference
figure
deltaTempC = tempC - avg24hTempC;
deltaTempC = reshape(deltaTempC, 24, 31).';

plot(1:24, mean(deltaTempC)), axis tight;
title('Mean temperature differential from 24 hour average');
xlabel('Hour of day (since midnight)');
ylabel('Temperature difference (\circC)');

%% Weighted Moving Average Filter
h = [1/2 1/2];  %for small n
binomialCoeff = conv(h,h);
for n = 1:4
    binomialCoeff = conv(binomialCoeff,h);
end

figure
fDelay = (length(binomialCoeff)-1)/2;
binomialMA = filter(binomialCoeff, 1, tempC);
plot(days, tempC, 'b', ...
     days-fDelay/24,  binomialMA, 'g');
axis tight;
legend('Hourly Temp','Binomial Weighted Average','location','best');
ylabel('Temp (\circC)');
xlabel('Time elapsed from Jan 1, 2011 (days)');
title('Logan Airport Dry Bulb Temperature (source: NOAA)');

%% Gaussian expansion filter
% adjust an exponentially weighted moving average filter by an alpha parameter between zero and one. A higher value of alpha will have less smoothing.
clf; alpha = 0.45;
exponentialMA = filter(alpha, [1 alpha-1], tempC);
plot(days, tempC, 'b', ...
     days-fDelay/24,  binomialMA, 'g', ...
     days-1/24,       exponentialMA,'r');

axis tight;
legend('Hourly Temp', ...
       'Binomial Weighted Average', ...
       'Exponential Weighted Average','location','best');
ylabel('Temp (\circC)');
xlabel('Time elapsed from Jan 1, 2011 (days)');
title('Logan Airport Dry Bulb Temperature (source: NOAA)');

%% open loop voltage experiment
load openloop60hertz
fs = 1000;
t = (0:numel(openLoopVoltage)-1) / fs;
plot(t,openLoopVoltage);
ylabel('Voltage (V)');
xlabel('Time (s)');
title('Open-loop Voltage Measurement');
% uniformly weighted moving average filter
plot(t,sgolayfilt(openLoopVoltage,1,17));
ylabel('Voltage (V)');
xlabel('Time (s)');
title('Open-loop Voltage Measurement');
legend('Moving average filter operating at 58.82 Hz', ...
       'location','southeast');
   
   %% linkdata test
x = 1:20;
y = rand(20,3);
area(x,y)
linkdata on
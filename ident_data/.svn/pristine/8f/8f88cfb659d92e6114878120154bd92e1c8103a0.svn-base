%% Estimating Simple Models from Real Laboratory Process Data
% This example shows how to develop and analyze simple models from a real
% laboratory process data. We start with a small description of the
% process, learn how to import the data to the toolbox and
% preprocess/condition it and then proceed systematically to estimate
% parametric and nonparametric models. Once the models have been identified
% we compare the estimated models and also validate the model to the actual
% output data from the experiment.

% Copyright 1986-2015 The MathWorks, Inc.

%% System Description
% This case study concerns data collected from a laboratory scale
% "hairdryer". (Feedback's Process Trainer PT326; See also page 525 in
% Ljung, 1999). The process works as follows: Air is fanned through a tube
% and heated at the inlet. The air temperature is measured by a
% thermocouple at the outlet. The input is the voltage over the heating
% device, which is just a mesh of resistor wires. The output is the outlet
% air temperature represented by the measured thermocouple voltage.

%% Setting up Data for Analysis
% First we load the input-output data to the MATLAB(R) Workspace.
load dryer2;

%%
% Vector |y2|, the output, contains 1000 measurements of the themocouple
% voltage which is proportional to the temperature in the outlet airstream.
% Vector |u2| contains 1000 input data points consisting of the voltage
% applied to the heater. The input was generated as a binary random
% sequence that switches from one level to the other with probability 0.2.
% The sample time is 0.08 seconds.
%
% The next step is to set up the data as an iddata object
dry = iddata(y2,u2,0.08);

%%
% To get information about the data, just type the name of the |iddata|
% object at the MATLAB command window:
dry

%%
% To inspect the properties of the above iddata object, use the |get|
% command:
get(dry)

%%
% For better book-keeping, it is good practice to give names to the input
% and output channels and Time units. These names would be propagated
% throughout the analysis of this iddata object:
dry.InputName = 'Heater Voltage';
dry.OutputName = 'Thermocouple Voltage';
dry.TimeUnit = 'seconds';
dry.InputUnit = 'V';
dry.OutputUnit = 'V';

%%
% Now that we have the data set ready, we choose the first 300 data points
% for model estimation.
ze = dry(1:300)

%% Preprocessing the Data
% Plot the interval from sample 200 to 300:
%
plot(ze(200:300));

%%
% *Figure 1:* A snapshot of the measured hair-dryer data.

%%
% From the above plot, it can be observed that the data is not zero mean.
% So let us remove the constant levels and make the data zero mean.
ze = detrend(ze);

%%
% The same data set after it has been detrended:
plot(ze(200:300)) %show samples from 200 to 300 of detrended data

%%
% *Figure 2:* Detrended estimation data.

%% Estimating Nonparametric and Parametric Models
% Now that the dataset has been detrended and there are no obvious
% outliers, let us first estimate the impulse response of the system
% by correlation analysis to get some idea of time constants and the
% like:

clf
mi = impulseest(ze); % non-parametric (FIR) model
showConfidence(impulseplot(mi),3); %impulse response with 3 standard 
                                   %deviations confidence region
%%
% *Figure 3:* Impulse response of the FIR model estimated using |ze|.

%%
% The shaded region marks a 99.7% confidence interval. There is a time
% delay (dead-time) of 3 samples before the output responds to the input
% (significant output outside the confidence interval). 

%%
% The simplest way to get started on a parametric estimation routine is to
% build a state-space model where the model-order is automatically
% determined, using a prediction error method. Let us estimate a model
% using the |ssest| estimation technique:
m1 = ssest(ze);

%%
% |m1| is a continuous-time identified state-space model, represented by an
% |idss| object. The estimation algorithm chooses 3 as the optimal order of
% the model. To inspect the properties of the estimated model, just enter
% the model name at the command window:
m1

%%
% The display suggests that the model is free-form (all entries of A, B and
% C matrices were treated as free parameters) and that the estimated model
% fits the data pretty well (over 90% fit). To retrieve the properties of
% this model, for example to obtain the |A| matrix of the discrete
% state-space object generated above, we can use the dot operator:
%
%       A = m1.a;
%
% See the "Data and Model Objects in System Identification Toolbox" example
% for more information regarding model objects. To find out which
% properties of the model object can be retrieved, use |get| command:
get(m1)

%%
% To fetch the values of the state-space matrices and their 1 std
% uncertainties, use the |idssdata| command:
[A,B,C,D,K,~,dA,dB,dC,dD,dK] = idssdata(m1)

%%
% The uncertainties are quite large even though the model fit the
% estimation data well. This is because the model is over-parameterized,
% that is, it has more free parameters than what could be uniquely
% identified. The variance of parameters in such cases is not well defined.
% However this does not imply that the model is unreliable. We can plot the
% time- and frequency-response of this plot and view the variance as
% confidence regions as discussed next. 

%% Analyzing the Estimated Model
%
% The Bode plot of the generated model can be obtained using the
% |bode| function as shown below:
h = bodeplot(m1);
%%
% *Figure 4:* Bode plot of estimated model.

%%
% Right-click on the plot and pick Characteristics->Confidence Region.
% Or, use the |showConfidence| command to view the variance of the
% response.
showConfidence(h,3) % 3 std (99.7%) confidence region 
%%
% *Figure 5:* Bode plot with 3 std confidence region.

%%
% Similarly, we can generate the step plot and its associated 3 std
% confidence region. We can compare the responses and associated variances
% of the parameteric model |m1| with that of the nonparametric model |mi|:
showConfidence(stepplot(m1,'b',mi,'r',3),3) 
%%
% *Figure 6:* Step plot of models |m1| and |mi| with confidence regions.

%%
% We can also consider the Nyquist plot, and mark uncertainty
% regions at certain frequencies with ellipses, corresponding to 3 
% standard deviations:
Opt = nyquistoptions;
Opt.ShowFullContour = 'off';
showConfidence(nyquistplot(m1,Opt),3)
%%
% *Figure 7:* Nyquist plot of estimated model showing the uncertainty regions at certain frequencies.

%%
% The response plots show that the estimated model |m1| is quite reliable. 

%% Estimating Models with a Prescribed Structure
% System Identification Toolbox can also be used to obtain a model with
% a prescribed structure. For example, a difference equation model with 2
% poles, 1 zero and 3 sample delays can be obtained using the |arx|
% function as shown below:
m2 = arx(ze,[2 2 3]);

%%
% To look at the model, enter the model name at the command window.
m2

%%
% A continuous time transfer function with 2 poles, one zero and 0.2 second
% transport delay can be estimated using the |tfest| command:
m3 = tfest(ze, 2, 1, 0.2)

%% Validating the Estimated Model to Experimental Output
% How good is an estimated model? One way to find out is to simulate it and
% compare the model output with measured output. Select a portion of the
% original data that was not used in building the model, say from samples
% 800 to 900. Once the validation data has been preprocessed, we use the
% |compare| function as shown below to view the quality of prediction:
zv = dry(800:900);   % select an independent data set for validation
zv = detrend(zv);    % preprocess the validation data
set(gcf,'DefaultLegendLocation','best')
compare(zv,m1);      % perform comparison of simulated output

%%
% *Figure 8:* Model's simulated response vs. validation data output.

%%
% It can be observed here that the agreement is very good. The "Fit" value
% shown is calculated as:
%
% |Fit = 100*(1 - norm(yh - y)/norm(y-mean(y)))|
%
% where |y| is the measured output (=|zv.y|), and |yh| is the output of the
% model |m1|.

%% Comparing Estimated Models
% To compare the performance of the models that we have estimated, for
% example |m1|, |m2| and |m3| with the validation data |zv|, we can again
% use the |compare| command:
compare(zv,m1,'b',m2,'r',m3,'c');

%%
% *Figure 9:* Comparing the responses of models |m1|, |m2|, |m3| on
% validation data set |ze|.

%%
% The pole-zero plots for the models can be obtained using |iopzplot|:

h = iopzplot(m1,'b',m2,'r',m3,'c');

%%
% *Figure 10:* Poles and zeros of the models |m1|, |m2| and |m3|.

%%
% The uncertainties in the poles and zeroes can also be obtained. In the
% following statement, '3' refers to the number of standard deviations.

showConfidence(h,3);

%%
% *Figure 11:* Pole-zero map with uncertainty regions.

%%
% The frequency functions above that are obtained from the models can be
% compared with one that is obtained using a non-parametric spectral
% analysis method (|spa|):
gs = spa(ze);

%%
% The |spa| command produces an IDFRD model. The |bode| function can again be
% used for a comparison with the transfer functions of the models obtained.
%
w = linspace(0.4,pi/m2.Ts,200);
opt = bodeoptions; opt.PhaseMatching = 'on';
bodeplot(m1,'b',m2,'r',m3,'c',gs,'g',w,opt);
legend('m1','m2','m3','gs')

%%
% *Figure 12:* Bode responses of |m1|, |m2| and |m3| compared against the
% non-parametric spectral estimation model |gs|.

%%
% The frequency responses from the three models/methods are very close.
% This indicates that this response is reliable.
%
% Also, a Nyquist plot can be analyzed with the uncertainty regions marked
% at certain frequencies:

showConfidence(nyquistplot(m1,'b',m2,'r',m3,'c',gs,'g'),3)
%%
% *Figure 13:* Nyquist plots of models |m1|, |m2|, |m3| and |gs|.

%%
% The non-parametric model |gs| exhibits the most uncertainty in response. 

%% Additional Information
% For more information on identification of dynamic systems with System
% Identification Toolbox visit the
% <http://www.mathworks.com/products/sysid/ System Identification Toolbox> product
% information page.

displayEndOfDemoMessage(mfilename)
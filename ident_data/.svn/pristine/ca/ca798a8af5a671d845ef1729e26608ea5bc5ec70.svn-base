1. Installation

 1.1  Create a directory, e.g., named as IDENT_DEMO_ACSE, in your local disk;

 1.2  Use one of the following approaches to copy the files to your destination directory created;
  (a) Use WinZip to extract the compressed file IDENT_DEMO.ZIP into your aim directory;
  (b) Use WinRAR to extract the compressed file IDENT_DEMO.RAR into your aim directory;
  (c) Copy directly all the files in IDENT_DEMO into your aim directory.

 1.3  Start MATLAB V6.x, and run a related file of this demo.

2. Objective
      The objective of this demonstration is to illustrate linear system identification
      techniques using MATLAB. Three MATLAB demonstrations are included:

 2.1  The identification of a gas furnace using
  (a) Correlation  analysis
  (b) Spectral analysis
  (c) Parameter estimation

 2.2  An example showing the effect of bias in parameter estimation.

 2.3  Two examples showing the time series prediction.

3. Brief discriptions on IDENT_DEMO.

 3.1 Data Files Used in Identification and Time Series Modelling 

 (1) GASDATA.DAT 
     The input and output data for a gas furnace system. There are two columns: the first
     column is the output of the system and the second one is the input. The data length
     is 296.
 
 (2) SUNSPOT_YEARLY.DAT
     The Wolf sunspot time series, annually recorded, two columns and the data length is 301.
     The first column is the year index (from 1700 to 2000), and the second column is the 
     yearly sunspot index.

 (3) SUNSPOT_MONTHly.DAT  
     The Wolf sunspot time series, monthly recorded, three columns and the data length is 3031.
     The first column is the year index (from 1747 to 2003), the second column is the months
     in each year (from Jul 1747 to March 2003), and the third is the monthly sunspot index.

 3.2 Programmes Used in Identification and Time Series Modelling

 (1) GAS_IDENT.m 
     This programme is used for the gas furnace system analysis and identification.
  
 (2) HEAT_EXCHANGER.m
     This programme is used for fitting AR or ARMAX models for the heat exchanger model based 
     on the band-limited white noise driven input-output data.

 (3) FREQ_ANALUSIS.m
     Frequency analysis for models (3.1), (3.2b) and (3.3).

 (4) SUNSPOT_AR.m
     This is used for fitting an AR model for sunspot time series.

 (5) SUNSPOT_ARMA.m
     Used for fitting an ARMA model for sunspot time series.

 (6) SIGNAL_DETECT.m
     For detecting a signal buried in a noise using cross correlation.

 3.3 Local Functions Used for Linear System Identification and Signal Processing
  
 (1) AR_ACSE.m
     Compute the parameters of an AR model whose order is given by a user.

 (2) ARX_ACSE.m
     Compute the parameters of an ARX model whose orders are given by a user.

 (3) ARMAX_ACSE.m
     Compute the parameters for an ARMA/ARMAX model whose order is given by a user.

 (4) CROSSCORR_ACSE.m
     Cross covariance function estimates.

 (5) CORELAT_ACSE.m
     Cross correlation function estimates.
  
 (6) PE_ACSE.m
     Compute the prediction errors.
  
 (7) RESID_ACSE.m
     Compute and test the residuals associated with a model.
  
 (8) SPECTRUM_ACSE.m
     Cross spectral estimates.
  
 (9) THETA_ACSE.m
     Specification for linear models and associated parameters. 
  
 (10)TH2POLY_ACSE.m
     Arrange the parameters of a given model as polynomial form.
  
 (11)XCORR_ACSE.m
     Normalized cross correlation function estimates.

    
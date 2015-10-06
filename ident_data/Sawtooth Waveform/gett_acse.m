function T=gett_acse(th);
%GETT_ACSE   Gets the sampling interval for a model.
%   T = GETT_ACSE(TH)
%
%   T: The sampling interval
%   TH: The model, defined in the THETA-format (See also THETA.)
%   See also SETT.

%   Copyright (c) 1986-98 by The MathWorks, Inc.
%   $Revision: 2.3 $  $Date: 1997/12/02 03:40:03 $

if nargin < 1
   disp('Usage: T = GETT(TH)')
   return
end

T=th(1,2);
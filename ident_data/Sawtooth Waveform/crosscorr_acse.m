function xycov = autocorr_acse(x,y,option)
%       Cross-covariance function estimates.
%	The cross-covariance is the cross-correlation function of
%	two sequences with their means removed.
%	COVARIAN(A,B), where A and B are length M vectors, returns the
%	length 2*M-1 cross-covariance sequence in a column vector.
%	COVARIAN(A), when A is a vector, is the auto-covariance sequence.
%	COVARIAN(A), when A is an M-by-N matrix, is a large matrix with
%	2*M-1 rows whose N^2 columns contain the cross-covariance
%	sequences for all combinations of the columns of A.
%	The zeroth lag of the output covariance is in the middle of the 
%	sequence, at element or row M.
%	By default, COVARIAN computes a raw covariance with no normalization.
%	COVARIAN(A,'biased') or COVARIAN(A,B,'biased) returns the "biased"
%	estimate of the cross-covariance function.  The biased estimate
%	scales the raw cross-covariance by 1/M.
%	COVARIAN(...,'unbiased') returns the "unbiased" estimate of the
%	cross-covariance function.  The unbiased estimate scales the raw
%	covariance by 1/(M-abs(k)), where k is the index into the result.
%	COVARIAN(...,'coeff') normalizes the sequence so that the
%	covariances at zero lag are identically 1.0.
%
%	References:
%	  [1] J.S. Bendat and A.G. Piersol, "Random Data:
%	      Analysis and Measurement Procedures", John Wiley
%	      and Sons, 1971, p.332.
%	  [2] A.V. Oppenheim and R.W. Schafer, Digital Signal 
%	      Processing, Prentice-Hall, 1975, pg 539.

[mx,nx] = size(x);
if nargin == 1
	xycov = correlat_acse(x-ones(mx,1)*mean(x));
elseif nargin == 2
	if isstr(y)
		xycov = correlat_acse(x-ones(mx,1)*mean(x),y);
	else
		[my,ny] = size(y);
		xycov = correlat_acse(x-ones(mx,1)*mean(x),y-ones(my,1)*mean(y));
	end
else
	[my,ny] = size(y);
	xycov = correlat_acse(x-ones(mx,1)*mean(x),y-ones(my,1)*mean(y),option);
end

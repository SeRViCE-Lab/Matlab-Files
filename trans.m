function y = trans(x, B, F)   %#codegen
% Preallocate output
x = linspace(600, 900, 10000);
y = zeros(size(x));

% Compute the transient on
y(1:(F+1)/2-1,:) = flipud(B((F-1)/2+2:end,:))*flipud(x(1:F,:));

% Compute the steady state output
ytemp = filter(B((F-1)./2+1,:),1,x);
y((F+1)/2:end-(F+1)/2+1,:) = ytemp(F:end,:);

% Compute the transient off
y(end-(F+1)/2+2:end,:) = flipud(B(1:(F-1)/2,:))*flipud(x(end-(F-1):end,:));

% Convert Y to the original shape of X
if isempty(DIM)
	y = shiftdim(y, -nshifts);
else
	y = ipermute(y,perm);
end
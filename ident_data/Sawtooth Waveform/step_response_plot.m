clc;
clear;
load step_response_31;
load step_response_32b;
load step_response_33;
x31=step_response_31(2,:);
x32b=step_response_32b(2,:);
x33=step_response_33(2,:);
nx31=length(x31);
nx32b=length(x32b);
nx33=length(x33);
nx=min([nx31,nx32b,nx33])-1;
plot(1:nx,x31(1:nx),1:nx,x32b(1:nx),'--'),
title('Unit Step Response of the Model (3.1)[solid line] and (3.2b)[dashed line]'), 
grid,
pause,
plot(1:nx,x31(1:nx),1:nx,x33(1:nx),'--'),
title('Unit Step Response of the Model (3.1)[solid line] and (3.3)[dashed line]'), 
grid,
pause,
close,


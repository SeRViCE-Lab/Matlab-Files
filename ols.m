%Olalekan Ogunmolu, Oklahoma City. Thanksgiving 2015.
%NARMAX Structure Identification using Orthogonal LS Estimators
clc; clear;
%cd('/home/lex/Documents/Matlab_Files/ident_data/UWGN')

data    = load('/home/lex/Documents/Matlab_Files/ident_data/UWGN/uwgn_data.mat');
data    = data.uwgn_data;

input   = data.input_ugwn;          %load some fictitious data
output  = data.output_ugwn;

%% Form Regressor matrix with na = 4 and nb = 4 as an example
% fit a linear-in-the parameters representation
N = length(input);

x = input;
y= output;
%parpool('local', 4);

%options = optimset(options, 'UseParallel', 'always');

for i = 1:N
    for j = 1:4
        for k = 5:8            
%          if(i == j)
%              Q(i,j) = -y(2);
%          else
            Q(i,j) = -y(j + 1);
%          end
%          if(i == k)
%              Q(i,k) = x(2);
%          else             
            Q(i,k) = x(k + 1);
%          end
        end
    end
end         

P = vertcat(output, input);
[W, A] = lu(P);                 %W = lower triangular, A = upper triangular

M =  length(P);

for k = 1:M
    for r = 1:(M-1)
        for m = 1:N
            num(m,:)    =   P(m,:) * W(r,:);
            Num(m,:)    =   sum(num);
            den         =   pow2 (W(r,:) );
            Den(m,:)    =   sum(den);            
        end
    end
end
%% compute g from the property of orthogonaility
Alpha = W' * W;

g = pinv(Alpha)  * W' * output;




clear all; clc;
x = [2,4,5,3,6]';

y = [6,2,0,1,7]';

u = x/norm(x);

v = y - (u'* y * u);

v = v/norm(v);

cost = x' * y/( norm(x)/norm(y) );

sint = sqrt( 1 - (cost)^2 );

R = eye(length(x))- u * u'- v * v' + [u v] * ...
      [cost -sint; sint cost] *[u v]';

fprintf('testing\n');

RRT = R*R'

R*x
function v = isodd(x)
%Tells you whether a number is odd or not.
%
% function v = isodd(x)
%
% v = (x/2 ~= round(x/2)) & (x == round(x));
v = (x/2 ~= round(x/2)) & (x == round(x));
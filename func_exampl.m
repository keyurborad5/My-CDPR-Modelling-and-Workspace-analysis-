clc
clear all
r=[2 5 8 66 77 41 2 9 54 66 25 33 225 11];
[a b c]=is_even(r);
[a; b; c]

function [y s w] = is_even(x)
y = ~mod(x,2);
s = ~mod(x,3);
w = ~mod(x,4);
%out=[y;s;w];
end

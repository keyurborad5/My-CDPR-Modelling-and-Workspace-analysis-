clc
clear all
A=[-5 5 -1;
    6 0 9];
for i=1:2
if any(A(i,:)<0) && any(A(i,:)>0)
    disp("YES")
end
end

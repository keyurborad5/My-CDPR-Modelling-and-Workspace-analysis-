clc
clear all
%% coordinates of first line segment
p1=[1,7,0.5];
q1=[4,5,0];
%% coordinates second line segment
p2=[2,5,-0.5];
q2=[5,6,0];
%% directional vectors of the segments
d1=-p1+q1;
d2=-p2+q2;
r=p1-p2;
a=dot(d1,d1);
b=dot(d1,d2);
e=dot(d2,d2);
c=dot(d1,r);
f=dot(d2,r);
denom=a*e-b*b;
if denom ~= 0
    s = clamp((b*f-c*e)/denom)
else s=0;
end

tnom=b*s+f;
if tnom<0
    t=0;
    s=clamp(-c/a)
elseif tnom>e
    t=1
    s=clamp((b-c)/a)
else 
    t=tnom/e;
    s=clamp((b*t-c)/a)
end
c1=p1+s*d1;
c2=p2+t*d2;
figure(1)
plot3([p1(1) q1(1)],[p1(2) q1(2)],[p1(3) q1(3)],[p2(1) q2(1)]...
    ,[p2(2) q2(2)],[p2(3) q2(3)],[c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)])
axis([0,9,0,9,-1,1])
grid on
dist = (dot(c1-c2,c1-c2))^0.5
function x = clamp(n)
if n<0
    x=0;
elseif n>1
    x=1;
else x=n;
end
end
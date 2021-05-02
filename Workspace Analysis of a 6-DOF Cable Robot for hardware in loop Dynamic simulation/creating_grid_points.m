clc
clear all
px=0:0.5:1;
py=0:0.5:1;
pz=0:0.5:1;
[X,Y,Z]=meshgrid(px,py,pz);
r=0;
for i=1:length(px)
    for j=1:length(py)
        for k=1:length(pz)
            P=[X(i,j,k);Y(i,j,k);Z(i,j,k)]
            r=r+1;
            plot3(P(1),P(2),P(3),".")
            hold on
            grid on
        end
    end
end
r
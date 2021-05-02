clc
clear all
%% First of all we will create a outer frame of CDPR 
% Frame of dimension (1mx1mx1m) and its corners are named with suffix A
% supplied coordinates are in world/Base frame
syms A;
A = [0,0,0;
    1,0,0;
    1,1,0;
    0,1,0;
    0,0,1;
    1,0,1;
    1,1,1;
    0,1,1];
%% Secondly will create an Object (cube) which we want to manipulate of dimension (0.1m,0.1m,0.1m)
%supplied coordinate will be in objects frame
%object frame is located at the COM of the object
syms B;
% B=0.01*[-5,5,-5;        %4
%         5,-5,5;         %6
%         -5,-5,-5;       %1
%         5,5,5;          %7
%          5,-5,-5;       %2
%          -5,5,5;        %8
%          5,5,-5;        %3
%          -5,-5,5];      %5
     %%%%%%%%%%%%%%%%%%%%%%%%%As given in papers
     B=0.01*[-10,5,-5;        %4
        10,-5,5;         %6
        -10,-5,-5;       %1
        10,5,5;          %7
         10,-5,-5;       %2
         -10,5,5;        %8
         10,5,-5;        %3
         -10,-5,5];      %5
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B = 0.01*[-5,-5,-5;
%     5,-5,-5;
%     5,5,-5;
%     -5,5,-5;
%     -5,-5,5;
%     5,-5,5;
%     5,5,5;
%     -5,5,5];
    
%% Now position of COM and Rotation matrix wrt base frame
syms R P phi psi theta;
Px=0.1:0.075:0.9;
% Px=0.5;
% Py=0.51;
%Pz=0.5;
Py=0.1:0.075:0.9;
Pz=0.1:0.075:0.9;
loop_no=0;
count=0;
[X,Y,Z]=meshgrid(Px,Py,Pz);
for i_1=1:length(Px)
    for j_1=1:length(Py)
        for k_1=1:length(Pz)
            P=[X(i_1,j_1,k_1);Y(i_1,j_1,k_1);Z(i_1,j_1,k_1)];
            loop_no=loop_no+1
% P = [0.5;0.5;0.1]; %position of COM of object in base frame
%fixed angle rotation x(psi),y(theta),z(phi)
%RXYZ = Rz(?)*Ry(?)*Rx(?)
psi=0*pi/180;phi=0*pi/180;theta=0;
R = [cos(phi)*cos(theta) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
    sin(phi)*cos(theta) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
    -sin(theta) cos(theta)*sin(psi) cos(theta)*cos(psi)] ;%Rotation of frame of Object wrt Base frame
%% Now vector representing the direction of the string and length of string is givrn as below
L = zeros(3,8);
Lm = zeros(8,1);
% a=[A(5,:);
%     A(5,:);
%     A(6,:);
%     A(6,:);
%     A(7,:);
%     A(7,:);
%     A(8,:);
%     A(8,:)];
a=[0,0.05,1;
    0.05,0,1;
    0.95,0,1;
    1,0.05,1;
    1,0.95,1;
    0.95,1,1;
    0.05,1,1;
    0,0.95,1];
for i=1:length(L)
    L(:,i)=a(i,:)'-(P+R*B(i,:)');
    Lm(i)=norm(a(i,:)'-(P+R*B(i,:)'));
end

L';
Lm;
%% Now to plot the cable from object end to frame end
O = zeros(3,8);
for i=1:length(O)
    O(:,i) = (P+R*B(i,:)');
end
% plot3([O(1,3) O(1,5)],[O(2,3) O(2,5)],[O(3,3) O(3,5)], ...
%     [O(1,5) O(1,7)],[O(2,5) O(2,7)],[O(3,5) O(3,7)], ...
%     [O(1,7) O(1,1)],[O(2,7) O(2,1)],[O(3,7) O(3,1)], ...
%     [O(1,1) O(1,3)],[O(2,1) O(2,3)],[O(3,1) O(3,3)], ...
%     [O(1,3) O(1,8)],[O(2,3) O(2,8)],[O(3,3) O(3,8)], ...
%     [O(1,8) O(1,2)],[O(2,8) O(2,2)],[O(3,8) O(3,2)], ...
%     [O(1,2) O(1,4)],[O(2,2) O(2,4)],[O(3,2) O(3,4)], ...
%     [O(1,4) O(1,6)],[O(2,4) O(2,6)],[O(3,4) O(3,6)], ...
%     [O(1,6) O(1,8)],[O(2,6) O(2,8)],[O(3,6) O(3,8)], ...
%     [O(1,6) O(1,1)],[O(2,6) O(2,1)],[O(3,6) O(3,1)], ...
%     [O(1,4) O(1,7)],[O(2,4) O(2,7)],[O(3,4) O(3,7)], ...
%     [O(1,2) O(1,5)],[O(2,2) O(2,5)],[O(3,2) O(3,5)], ...
%     'LineWidth',2)
%     hold on
% plot3(P(1),P(2),P(3),'o','MarkerFaceColor','red','MarkerSize',10)
% axis on
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% for i=1:8
%     plot3([a(i,1) O(1,i)],[a(i,2) O(2,i)],[a(i,3) O(3,i)],'b','LineWidth',1)
%     hold on
%     grid on
%     axis([0 1 -0 1 -0 1])
% end
% 
% hold off 
%% Static and kinematic model
% Wrench matrics
unit_v=zeros(8,3);
for i=1:8
unit_v(i,:)=L(:,i)'/norm(L(:,i));
end

for i=1:8
% cross_prod(:,i) = cross(O(:,i),unit_v(i,:)');
cross_prod(:,i) = cross(R*B(i,:)',unit_v(i,:)');
end

W = [unit_v';cross_prod];
r_k = rank(W);
% W=[-0.5946    0.5946   -0.5946    0.5946    0.5946   -0.5946    0.5946   -0.5946;
%     0.4460   -0.4460   -0.4460    0.4460   -0.4460    0.4460    0.4460   -0.4460;
%    -0.6690    0.6690   -0.6690    0.6690   -0.6690    0.6690   -0.6690    0.6690;
%    -0.6690    0.4460         0    0.2230         0    0.2230   -0.6690    0.4460;
%    0   -0.0743   0   -0.0743    0.6690   -0.5946    0.6690   -0.5946;
%     0.5946   -0.4460    0   -0.1487   -0.4460    0.5946   -0.1487         0];
Arr=1:8;
index =nchoosek(Arr,5);
for i=1:length(index)
W_65(:,:,i)=W(:,index(i,:));
c_k(:,i)=null(W_65(:,:,i)');
for j=1:6
if abs(c_k(j,i))<=1e-3
    c_k(j,i)=0;
end
end 
end
c_l=-c_k;
c=[c_k c_l];
Tmax=100;
Tmin=2;
for i=1:length(c_k)
    I_plus=0;
    I_minus=0;
    for j=1:length(W)
        if c_k(:,i)'*W(:,j)>0
            I_plus=I_plus+c_k(:,i)'*W(:,j);
        else 
            I_minus=I_minus+c_k(:,i)'*W(:,j);
        end
        d_k(i)=Tmax*I_plus+Tmin*I_minus;
        d_l(i)=-Tmax*I_minus-Tmin*I_plus;
    end
end
d=[d_k d_l];
Mmax=2;
g=9.81;
for j=1:length(d)
    r1(j)=(d(j)/(Mmax*g)-c(3,j))/((c(4,j)^2+c(5,j)^2)^0.5);
end
r1_min=min(r1);
if r1_min>=0.05
     figure(1)
    plot3(P(1),P(2),P(3),".r")
    axis([-0.1 1.1 -0.1 1.1 -0.1 1.1])
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('z-axis');
    grid on
    hold on
else count=count+1;
end
        end
    end
end
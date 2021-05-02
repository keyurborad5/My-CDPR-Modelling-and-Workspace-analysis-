clc
clear all
%% First of all we will create a outer frame of CDPR 
% Frame of dimension (1mx1mx1m) and its corners are named with suffix A
% supplied coordinates are in world/Base frame
tic
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
     B=0.01*[-10,20,-5;        %4
        10,-20,5;         %6
        -10,-20,-5;       %1
        10,20,5;          %7
         10,-20,-5;       %2
         -10,20,5;        %8
         10,20,-5;        %3
         -10,-20,5];      %5
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
Px=0.05:0.1:0.9;
Py=0.05:0.1:0.9;
Pz=0.05:0.1:0.9;
[X,Y,Z]=meshgrid(Px,Py,Pz);
for i_1=1:length(Px)
    for j_1=1:length(Py)
        for k_1=1:length(Pz)
            P=[X(i_1,j_1,k_1);Y(i_1,j_1,k_1);Z(i_1,j_1,k_1)];
%P = [0.5;0.7;0.5]; %position of COM of object in base frame
%fixed angle rotation x(psi),y(theta),z(phi)
%RXYZ = Rz(phi)*Ry(theta)*Rx(psi)
psi=0*pi/180;phi=-15*pi/180;theta=0*pi/180; %provide angles in radians
R = [cos(phi)*cos(theta) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
    sin(phi)*cos(theta) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
    -sin(theta) cos(theta)*sin(psi) cos(theta)*cos(psi)]; %Rotation of frame of Object wrt Base frame
%% Now vector representing the direction of the string and length of string is givrn as below
L = zeros(3,8);
Lm = zeros(8,1);
% a=[A(1,:);
%    A(5,:);
%    A(2,:);
%    A(6,:);
%     A(3,:);
%     A(7,:);
%     A(4,:);
%     A(8,:)];
%**************************
a=[A(4,:);
    A(6,:);
    A(1,:);
    A(7,:);
    A(2,:);
    A(8,:);
    A(3,:);
    A(5,:)];

%************************
for i=1:length(L)
    L(:,i)=a(i,:)'-(P+R*B(i,:)');
    Lm(i)=norm(a(i,:)'-(P+R*B(i,:)'));
end
L';
Lm;
%% Central platform vertices in global frame
O = zeros(3,8);
for i=1:length(O)
    O(:,i) = (P+R*B(i,:)');
end
%% Static and kinematic model
% Wrench matrics
unit_v=zeros(8,3);
for i=1:8
unit_v(i,:)=L(:,i)'/norm(L(:,i));
end
for i=1:8
% cross_prod(:,i) = cross(O(:,i),unit_v(i,:)'); this is incorect pls nnote
cross_prod(:,i) = cross(R*B(i,:)',unit_v(i,:)');
end
W = [unit_v';cross_prod];
r_k = rank(W);
for i=1:6
    for j=1:8
        if abs(W(i,j))<=1e-4
            W(i,j)=0;
        end
    end
end
%% Workspace Algorithmn
rk = rank(W);
Wnew=W(:,1:6);
% invWn=inv(Wnew);
% for i=1:6
%     for j=1:6
%         if abs(invWn(i,j))<=1e-3
%             invWn(i,j)=0;
%         end
%     end
% end
A_78=Wnew\W(:,7:8);
for i=1:6
    for j=1:2
        if abs(A_78(i,j))<=1e-2
            A_78(i,j)=0;
        end
    end
end
%% For printing DOT if the workspace position is feasible
co_nt=0;
for i=1:6
    if A_78(i,1)<0 || A_78(i,2)<0
        co_nt=co_nt+1;
        if co_nt==6
            figure(1)
            plot3(P(1),P(2),P(3),".r")
            axis([-0.1 1.1 -0.1 1.1 -0.1 1.1])
            xlabel('x-axis');
            ylabel('y-axis');
            zlabel('z-axis');
            grid on
            hold on
        end
    end
end
end
end
end
toc
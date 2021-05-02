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
% B = 0.01*[-10,-20,-5;
%     10,-20,-5;
%     10,20,-5;
%     -10,20,-5;
%     -10,-20,5;
%     10,-20,5;
%     10,20,5;
%     -10,20,5];
    
%% Now position of COM and Rotation matrix wrt base frame
syms R P phi psi theta;
P = [0.5;0.5;0.1]; %position of COM of object in base frame
%fixed angle rotation x(psi),y(theta),z(phi)
%RXYZ = Rz(phi)*Ry(theta)*Rx(psi)
psi=45*pi/180;phi=0*pi/180;theta=0; %provide angles in radians
R = [cos(phi)*cos(theta) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
    sin(phi)*cos(theta) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
    -sin(theta) cos(theta)*sin(psi) cos(theta)*cos(psi)] %Rotation of frame of Object wrt Base frame
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
%*****************************
a=[0,0.05,1;
    0.05,0,1;
    0.95,0,1;
    1,0.05,1;
    1,0.95,1;
    0.95,1,1;
    0.05,1,1;
    0,0.95,1];
%***************************
%**************************
% a=[A(4,:);
%     A(6,:);
%     A(1,:);
%     A(7,:);
%     A(2,:);
%     A(8,:);
%     A(3,:);
%     A(5,:)];

%************************
for i=1:length(L)
    L(:,i)=a(i,:)'-(P+R*B(i,:)');
    Lm(i)=norm(a(i,:)'-(P+R*B(i,:)'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=[L ([P(1);P(2);0]-(P+R*[0;0;-0.05]))]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L';
Lm;
%% Now to plot the cable from object end to frame end
O = zeros(3,8);
for i=1:length(O)
    O(:,i) = (P+R*B(i,:)');
end
plot3([O(1,3) O(1,5)],[O(2,3) O(2,5)],[O(3,3) O(3,5)], ...
    [O(1,5) O(1,7)],[O(2,5) O(2,7)],[O(3,5) O(3,7)], ...
    [O(1,7) O(1,1)],[O(2,7) O(2,1)],[O(3,7) O(3,1)], ...
    [O(1,1) O(1,3)],[O(2,1) O(2,3)],[O(3,1) O(3,3)], ...
    [O(1,3) O(1,8)],[O(2,3) O(2,8)],[O(3,3) O(3,8)], ...
    [O(1,8) O(1,2)],[O(2,8) O(2,2)],[O(3,8) O(3,2)], ...
    [O(1,2) O(1,4)],[O(2,2) O(2,4)],[O(3,2) O(3,4)], ...
    [O(1,4) O(1,6)],[O(2,4) O(2,6)],[O(3,4) O(3,6)], ...
    [O(1,6) O(1,8)],[O(2,6) O(2,8)],[O(3,6) O(3,8)], ...
    [O(1,6) O(1,1)],[O(2,6) O(2,1)],[O(3,6) O(3,1)], ...
    [O(1,4) O(1,7)],[O(2,4) O(2,7)],[O(3,4) O(3,7)], ...
    [O(1,2) O(1,5)],[O(2,2) O(2,5)],[O(3,2) O(3,5)], ...
    'LineWidth',2)
    hold on
plot3(P(1),P(2),P(3),'o','MarkerFaceColor','red','MarkerSize',10)
axis on
xlabel('x')
ylabel('y')
zlabel('z')
hold on
for i=1:8
    plot3([a(i,1) O(1,i)],[a(i,2) O(2,i)],[a(i,3) O(3,i)],'b','LineWidth',1)
    hold on
    grid on
    axis([-0.5 1.5 -0.5 1.5 -0.3 1.5])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot3([P(1) P(1)],[P(2) P(2)],[P(3)-0.05 0],'g','LineWidth',1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off 
%% Static and kinematic model
% Wrench matrics
unit_v=zeros(8,3);
for i=1:8
unit_v(i,:)=L(:,i)'/norm(L(:,i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unit_v = [unit_v; L(:,9)'/norm(L(:,9))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:8
% cross_prod(:,i) = cross(O(:,i),unit_v(i,:)');
cross_prod(:,i) = cross(R*B(i,:)',unit_v(i,:)');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cross_prod = [cross_prod cross([P(1)-P(1);P(2)-P(2);P(3)-0.05],unit_v(9,:)')]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = [unit_v';cross_prod]
r_k = rank(W)
%% Wrench applied to the platform due to gravity Wg
m_p = 1;%mass of platform as 1KG
MS_p = R*m_p*P;
g=[0;0;-9.81];
Wg = [m_p*eye(3);
        [0 -MS_p(3) MS_p(2);
        MS_p(3) 0 -MS_p(1);
        -MS_p(2) MS_p(1) 0]]*g;
%% Finding tension in the string
% W*t+We+Wg=0
% here t=pinv(W)*(-Wg)Moore-Penrose pseudoinverse
t = pinv(W)*(-Wg-[0;0;9;4.5;-4.5;0]); %+ null(W,'r')*[2;1]
norm(t);
%% Workspace Algorithmn
rk = rank(W)
Wnew=W(:,1:6);
invWn=inv(Wnew);
A=invWn*W
for i=1:6
    for j=1:9
        if abs(A(i,j))<=1e-3
            A(i,j)=0;
        end
    end
end
A
% E=zeros(6,1);
% for i=1:8
%     E=+W(:,i);
% end
% E
toc
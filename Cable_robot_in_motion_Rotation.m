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
B=0.01*[-5,5,-5;        %4
        5,-5,5;         %6
        -5,-5,-5;       %1
        5,5,5;          %7
         5,-5,-5;       %2
         -5,5,5;        %8
         5,5,-5;        %3
         -5,-5,5];      %5
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
% here trajectory is generated using Transform matrix for rotation and
% translation simultaneously
t0 = trvec2tform([0.2 0.5 0.2])*axang2tform([1 1 1 0]);
tF = trvec2tform([0.7 0.5 0.7])*axang2tform([0 0 1 -pi/3]);
tInterval = [0 1];
tvec = 0:0.01:1;
[tfInterp, v1, a1] = transformtraj(t0,tF,tInterval,tvec);
q_t  = tform2trvec(tfInterp);
q = q_t';
R =  tform2rotm(tfInterp);
% wpts = [0.2 0.3 0.4 0.5 0.6 0.7;
%         0.5 0.5 0.5 0.5 0.5 0.5;
%         0.5 0.5 0.5 0.5 0.5 0.5];
% tpts = [0 1 2 3 4 5];
% tvec = 0:0.1:5;
% [q, qd, qdd, pp] = quinticpolytraj(wpts, tpts, tvec); % Here we are using 5th order polynomial to fit into given waypoints(Here motion is only in 3 Axis translation)
%P = [0.5;0.5;0.3]; %position of COM of object in base frame
%fixed angle rotation x(psi),y(theta),z(phi)
%RXYZ = Rz(?)*Ry(?)*Rx(?)
% psi=0;phi=0;theta=0;
% R = [cos(phi)*cos(theta) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
%     sin(phi)*cos(theta) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
%     -sin(theta) cos(theta)*sin(psi) cos(theta)*cos(psi)] %Rotation of frame of Object wrt Base frame
%% Now vector representing the direction of the string and length of string is givrn as below
L = zeros(3,8,length(q));
Lm = zeros(8,length(q));
a=[A(5,:);
    A(5,:);
    A(6,:);
    A(6,:);
    A(7,:);
    A(7,:);
    A(8,:);
    A(8,:)];
for j= 1:length(q)
    for i=1:8
        L(:,i,j)=a(i,:)'-(q(:,j)+R(:,:,j)*B(i,:)');
        Lm(i,j)=norm(a(i,:)'-(q(:,j)+R(:,:,j)*B(i,:)'));
    end
end

Lm;
%% Now to plot the cable from object end to frame end
O = zeros(3,8,length(q));
for j=1:length(q)
    for i=1:8
        O(:,i,j) = (q(:,j)+R(:,:,j)*B(i,:)');
    end
    
    figure(1)
    plot3([O(1,3,j) O(1,5,j)],[O(2,3,j) O(2,5,j)],[O(3,3,j) O(3,5,j)], ...
    [O(1,5,j) O(1,7,j)],[O(2,5,j) O(2,7,j)],[O(3,5,j) O(3,7,j)], ...
    [O(1,7,j) O(1,1,j)],[O(2,7,j) O(2,1,j)],[O(3,7,j) O(3,1,j)], ...
    [O(1,1,j) O(1,3,j)],[O(2,1,j) O(2,3,j)],[O(3,1,j) O(3,3,j)], ...
    [O(1,3,j) O(1,8,j)],[O(2,3,j) O(2,8,j)],[O(3,3,j) O(3,8,j)], ...
    [O(1,8,j) O(1,2,j)],[O(2,8,j) O(2,2,j)],[O(3,8,j) O(3,2,j)], ...
    [O(1,2,j) O(1,4,j)],[O(2,2,j) O(2,4,j)],[O(3,2,j) O(3,4,j)], ...
    [O(1,4,j) O(1,6,j)],[O(2,4,j) O(2,6,j)],[O(3,4,j) O(3,6,j)], ...
    [O(1,6,j) O(1,8,j)],[O(2,6,j) O(2,8,j)],[O(3,6,j) O(3,8,j)], ...
    [O(1,6,j) O(1,1,j)],[O(2,6,j) O(2,1,j)],[O(3,6,j) O(3,1,j)], ...
    [O(1,4,j) O(1,7,j)],[O(2,4,j) O(2,7,j)],[O(3,4,j) O(3,7,j)], ...
    [O(1,2,j) O(1,5,j)],[O(2,2,j) O(2,5,j)],[O(3,2,j) O(3,5,j)], ...
    'LineWidth',2)
    hold on
    plot3(q(1,j),q(2,j),q(3,j),'o','MarkerFaceColor','red','MarkerSize',10)
    axis on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
    for i=1:8
        plot3([a(i,1) O(1,i,j)],[a(i,2) O(2,i,j)],[a(i,3) O(3,i,j)],'b','LineWidth',1)
        hold on
        grid on
        axis([0 1 0 1 -0.3 1.5])
    end
    hold off
end
%% Static and kinematic model
% Wrench matrics
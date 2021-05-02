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
    

%% Now vector representing the direction of the string and length of string is givrn as below

a=[A(5,:);
   A(5,:);
   A(6,:);
   A(6,:);
   A(7,:);
   A(7,:);
   A(8,:);
   A(8,:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a=[A(4,:);
%     A(6,:);
%     A(1,:);
%     A(7,:);
%     A(2,:);
%     A(8,:);
%     A(3,:);
%     A(5,:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now position of COM and Rotation matrix wrt base frame
syms R P phi psi theta;
% Px=0.5;
% Py=0.5;
% Pz=0.5;
Px=-0.1:0.1:1.1;
Py=-0.1:0.1:1.1;
Pz=0:0.1:1;
loop_no=0;
[X,Y,Z]=meshgrid(Px,Py,Pz);
for i_1=1:length(Px)
    for j_1=1:length(Py)
        for k_1=1:length(Pz)
            P=[X(i_1,j_1,k_1);Y(i_1,j_1,k_1);Z(i_1,j_1,k_1)];
            loop_no=loop_no+1
            
% P = [0.5;0.5;0.5]; %position of COM of object in base frame
%fixed angle rotation x(psi),y(theta),z(phi)
%RXYZ = Rz(?)*Ry(?)*Rx(?)
psi=0*pi/180;phi=0*pi/180;theta=0*pi/180; %provide angles in radians
R = [cos(phi)*cos(theta) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
    sin(phi)*cos(theta) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
    -sin(theta) cos(theta)*sin(psi) cos(theta)*cos(psi)]; %Rotation of frame of Object wrt Base frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = zeros(3,8);
Lm = zeros(8,1);
for i=1:length(L)
    L(:,i)=a(i,:)'-(P+R*B(i,:)');
    Lm(i)=norm(a(i,:)'-(P+R*B(i,:)'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=[L ([P(1);P(2);0]-(P+R*[0;0;-0.05]))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L';
Lm;
%% Now to plot the cable from object end to frame end
O = zeros(3,8);
for i=1:length(O)
    O(:,i) = (P+R*B(i,:)');
end
%% Static and kinematic model
% Wrench matrics
cross_prod = zeros(3,8);
unit_v=zeros(8,3);
for i=1:8
unit_v(i,:)=L(:,i)'/norm(L(:,i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unit_v = [unit_v; L(:,9)'/norm(L(:,9))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:8
cross_prod(:,i) = cross(R*B(i,:)',unit_v(i,:)');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cross_prod = [cross_prod cross([P-[P(1);P(2);0.05]],unit_v(9,:)')];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = [unit_v';cross_prod];
for i=1:6
    for j=1:9
        if abs(W(i,j))<=1e-4
            W(i,j)=0;
        end
    end
end
W;
r_k = rank(W);
%% Workspace determination Algorithm
Arr=1:9;
index = nchoosek(Arr,5);
for j = 1:length(index)        % Now we want to fill the remaing 3 columns that we excluded earlier 
    k=0;
    for i=1:9
        if index(j,:)~=Arr(1,i)
            k=k+1;
            d(j,k)=Arr(1,i);
        end
    end
end
com_index = [index d]; %here i have now 5 combi of index from 8 index and remaining are appended
k=0;
for i = 1:length(com_index)
    W_temp(:,:,i)= W(:,com_index(i,:));     %using index to rearrange columns accord.
    if rank(W_temp(:,1:5,i))==5
        k=k+1;
        Wnew(:,:,k)=W_temp(:,:,i);      %after checking rank of first5 cols W_temp is assigned to form Wnew
    end
end
for i =1:length(Wnew)
    W_65=Wnew(:,1:5,i);             %formed new mat. W_65 from Wnew (first 5 cols)
    ran_k(i,:) = rank(Wnew(:,1:5,i));
    for j=1:6
        W_65(j,:)=[];               %eleminated an row to get 5X5 mat. to compute det. for finding component of Nrml vector
        r_k(j,i)=rank(W_65);
        Nrml(j,i)=(-1)^(j+1)*det(W_65);
        W_65=Wnew(:,1:5,i);
    end
end
for i=1:6
    for j=1:length(Nrml)
        if abs(Nrml(i,j))<=1e-4
            Nrml(i,j)=0;
        end
    end
end
for i = 1:length(Nrml)
    s(i,:) = Nrml(:,i)'*Wnew(:,6:9,i); %carried out dot product of Nrml and remaining cols of corros. Wnew
end
f=0;
for i = 1:length(s)  % this loop is for setting low order value to zeros
%     for j = 1:3
%         if abs(s(i,j))<=0.0001
%             s(i,j) = 0;
%         end
%     end
    if any(s(i,:)>0) && any(s(i,:)<0)
        f=f+1;
    end
end
if f==length(s)
     figure(1)
    plot3(P(1),P(2),P(3),".r")
    axis([-0.15 1.15 -0.15 1.15 -0.1 1.1])
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('z-axis');
    grid on
    hold on
end
        end
    end
end
time_elapsed=toc
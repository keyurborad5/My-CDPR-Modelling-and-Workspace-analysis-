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
P = [0.2;0.2;0.6]; %position of COM of object in base frame
%fixed angle rotation x(psi),y(theta),z(phi)
%RXYZ = Rz(?)*Ry(?)*Rx(?)
psi=0;phi=10*pi/180;theta=0; %provide angles in radians
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
%    A(3,:);
%    A(7,:);
%    A(4,:);
%    A(8,:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=[A(4,:);
    A(6,:);
    A(1,:);
    A(7,:);
    A(2,:);
    A(8,:);
    A(3,:);
    A(5,:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
hold off 
%% Static and kinematic model
% Wrench matrics
unit_v=zeros(8,3);
for i=1:8
unit_v(i,:)=L(:,i)'/norm(L(:,i));
end
for i=1:8
% cross_prod(:,i) = cross(O(:,i),unit_v(i,:)');this is incorect pls note
cross_prod(:,i) = cross(R*B(i,:)',unit_v(i,:)');
end
W = [unit_v';cross_prod]
for i=1:6
    for j=1:8
        if abs(W(i,j))<=1e-3
            W(i,j)=0;
        end
    end
end
W
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
% clc
% clear all
% W = [-0.5637   -0.5851    0.5108    0.6353    0.5637    0.5851   -0.5108   -0.6353;
%    -0.5108   -0.6353   -0.5637   -0.5851    0.5108    0.6353    0.5637    0.5851;
%     0.6490    0.5041    0.6490    0.5041    0.6490    0.5041    0.6490    0.5041;
%     0.5108    0.6353    0.5637    0.5851    0.1382   -0.1312    0.0853   -0.0810;
%    -0.5637   -0.5851   -0.1382    0.1312   -0.0853    0.0810   -0.5108   -0.6353;
%          0         0   -0.5637   -0.5851   -0.0529    0.0502    0.5108    0.6353];
Arr=1:8;
index = nchoosek(Arr,5);
for j = 1:length(index)        % Now we want to fill the remaing 3 columns that we excluded earlier 
    k=0;
    for i=1:8
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
        if abs(Nrml(i,j))<=1e-3
            Nrml(i,j)=0;
        end
    end
end
for i = 1:length(Nrml)
    s(i,:) = Nrml(:,i)'*Wnew(:,6:8,i); %carried out dot product of Nrml and remaining cols of corros. Wnew
end
f=0;
for i = 1:length(s)  % this loop is for setting low order value to zeros
    for j = 1:3
        if abs(s(i,j))<=0.001
            s(i,j) = 0;
        end
    end
    if s(i,1)~=0 && s(i,2)~=0 && s(i,3)~=0
        f=f+1;
        s_new(f,:)=s(i,:);
    end
end
% f=0;
% for i = 1:length(s)  % for finding out rows which doesnot contain any zeros in them
% if s(i,1)~=0 & s(i,2)~=0 & s(i,3)~=0
% f=f+1;
% s_new(f,:)=s(i,:);
% end
% end
f=0;
for i = 1: length(s) % checking for sign change in rows with no zeros 
    if s(i,2)/s(i,1)>=0 && s(i,3)/s(i,1)>=0
        f=f+1;
        s_sign(f,:)=s(i,:);
    end
end
time_elapsed=toc
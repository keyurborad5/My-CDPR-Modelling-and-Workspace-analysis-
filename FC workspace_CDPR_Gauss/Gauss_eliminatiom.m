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
P = [0.7;0.5;0.5]; %position of COM of object in base frame
%fixed angle rotation x(psi),y(theta),z(phi)
%RXYZ = Rz(phi)*Ry(theta)*Rx(psi)
psi=0;phi=0;theta=0; %provide angles in radians
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
% cross_prod(:,i) = cross(O(:,i),unit_v(i,:)');this is incorrect plis note
cross_prod(:,i) = cross(R*B(i,:)',unit_v(i,:)');
end
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
%% Gauss elimination
%W = [0 1 1;2 0 3;1 1 1];
A_p=1:8;              %creating a row array with 8 elements 1 to 8
a_p = perms(A_p);      %permutation of nPr = 8P8 and storing all different arrangments(i.e 40320) in rows 
b= a_p(:,1:5);        %selecting first 5 columns of matix 'a'
c = unique(b,'rows'); %as only first 5 columns are selected there will be dublicates of the rows, hence elimating them. this is sililar to 8P5 and wrting all possible arrangment in rows
d=zeros(6720,3);        % created a matrix of zeros of size 6720X3 for sack of declaration to use ahead.
for j = 1:length(c)        % Now we want to fill the remaing 3 columns that we excluded earlier 
    k=0;
    for i=1:8
        if c(j,:)~=A_p(1,i)
            k=k+1;
            d(j,k)=A_p(1,i);
        end
    end
end
ind = [c d];    % hence made 8P5=6720 arrangments and appended remaining 3 columns
for i=1:6720
    W(:,:,i)=W(:,ind(i,:)); % created W[6 X 8] for all possible arrangment and found out above. W[6 X 8 X 6720]
end
[n,~]=size(W);      % from here we start the method of gauss elimination for in all 6720 matrix
f=0;
W78 = W(:,:,78);
for j=1:6720
    all_rank(j)=rank(W(:,:,j));
    if rank(W(:,:,j))==6    % we egnore matrix which are not full rank that is 6
        count = 0;
        for i=1:n-1
            
            if abs(W(i,i,j))<=1e-6  % this case runs when your pivot point turns to be zero 
                
                for q= i+1:n
                    if abs(W(q,i,j))>=1e-6      %we swap this row with the rows below it if they are not zero
                        count = count+1;
                        temp = W(i,:,j);
                        W(i,:,j)=W(q,:,j);
                        
                        W(q,:,j)=temp;
                        break
                    end
%               W(:,:,j)=[W(:,:,j);W(i,:,j)];
%               W(i,:,j)=[];
                end
            end
            if count~=0| abs(W(i,i,j))>=1e-6
                   W(i,:,j)=W(i,:,j)/W(i,i,j);  %actual code for Gauss elimination
                   m=W(i+1:n,i,j)/W(i,i,j);
                   W(i+1:n,:,j)=W(i+1:n,:,j)-m*W(i,:,j);
            end
        end
        f=f+1;
        oneD(f,:) = W(n,6:8,j); %this stores the last rows last 3 columns terms which we need to eval
    end
end
% For force closure to satisfy we need to satisy convex hull and CH is
% carried out with help of Gaussian Elimination. In this meathod the last 3
% terms of last row contains values of opposite sign then convex hull
% encloses the origin and hence satisfy Force closure and that pose is
% achiveable and valid in real life.
% here we have considered that if zero is present in the last 3 elements
% than also it encloses the origin in convex hull(bt its not mentioned in
% the r-paper)
for i = 1:6720  % this loop is for setting low order value to zeros
for j = 1:3
if abs(oneD(i,j))<=0.001
oneD(i,j) = 0;
end
end
end
f=0;
for i = 1:6720  % for finding out rows which doesnot contain any zeros in them
if oneD(i,1)~=0 & oneD(i,2)~=0 & oneD(i,3)~=0
f=f+1;
nz(f,:)=oneD(i,:);
end
end
f=0;
for i = 1: length(nz) % checking for sign change in rows with no zeros 
    if nz(i,2)/nz(i,1)>0 & nz(i,3)/nz(i,1)>0
        f=f+1;
        s_sign(f,:)=nz(i,:);
    end
end
timeElapsed = toc
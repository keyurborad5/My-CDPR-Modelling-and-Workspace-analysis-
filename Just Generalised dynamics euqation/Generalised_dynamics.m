clc
clear all
%% here need to mention first of all all the symbolic variables
syms th_1 th_2 th_3 a_1 a_2 a_3 m1 m2 m3 Ic1 Ic2 Ic3 real
mass=[m1 m2 m3];
Ic=[Ic1 Ic2 Ic3];
q=[th_1 th_2 th_3];
M=0;
%% INput of DH parameters (A Alpha D Theta)
prompt = "ENTER THE DH parameter (a alpha d theta) : ";
dh_para = input(prompt)
%% Input of Angular Jacobian matrix of the end effector
ang_j = "enter the end effector angular jacobian: ";
Jacobian_w = input(ang_j)
[n,~]=size(dh_para);
% for i=1:n
%     mass(i)=input ("enter the mass: ")
%     MOI(i)=input("Enter the MOI: ")
% end
%% Calling for function which converts the dh paramenters to Transformation matrix
for i=1:n
Tr_mat(:,:,i) = Homo(dh_para(i,:));
end
%% Multiplying matrix to get 0T2=0T1*1T2 etc.
for i=1:n-1
TR(:,:,i)= Tr_mat(:,:,i)*Tr_mat(:,:,i+1)
end
%% position of links
Pc_1=Tr_mat(1:3,4,1)
for i=1:n-1
    Pc(:,i)=TR(1:3,4,i)
end
Pc = [Pc_1 Pc];
Pc=simplify(Pc)

%% Jacobian
for i=1:n
    for j=1:n
        Jac_v(:,j,i)=[diff(Pc(:,i),q(j))]; %Linear velocity jacobian
    end
end
Jac_v=simplify(Jac_v)
for i=1:n
    Jac_w(:,:,i)= [Jacobian_w(:,1:i) zeros(3,n-i)] % this is arranging the Angular Velocity Jacobian
end
%% Getting mass matrix by adding all the terms with each other
for i=1:n
M = M + mass(i)*Jac_v(:,:,i)'*Jac_v(:,:,i)+Ic(i)*Jac_w(:,:,i)'*Jac_w(:,:,i);
end
M= simplify(M)
%% function of dh para to Homoeneous transformation
function TR = Homo(x)
TR = [cos(x(4)) -sin(x(4))*cos(x(2)) sin(x(4))*sin(x(2)) x(1)*cos(x(4));
    sin(x(4)) cos(x(4))*cos(x(2)) -sin(x(2))*cos(x(4)) x(1)*sin(x(4));
    0 sin(x(2)) cos(x(2)) x(3);
    0 0 0 1];
end

function k_ele=PlaneFrameElementStiffness(A,E,I,L,angle)
%PlaneFrameElementStiffness This function returns the element
%element with modules of elasticity E,
%cross-sectional area A, moment of
%inertia I, length L, and angle
%theta (in degrees).
%The size of the element stiffness
%matrix is 6 x 6.


c=cosd(angle);
s=sind(angle);
C1=A*E/L;
C2=E*I/L^3;
k_local=[C1,0,0,-C1,0,0;
         0,12*C2,6*C2*L,0,-12*C2,6*C2*L;
         0,6*C2*L,4*C2*L*L,0,-6*C2*L,2*C2*L*L;
         -C1,0,0,C1,0,0;
         0,-12*C2,-6*C2*L,0,12*C2,-6*C2*L;
         0,6*C2*L,2*C2*L*L,0,-6*C2*L,4*C2*L*L;]
T=[ c s 0 0 0 0;
   -s c 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 c s 0;
    0 0 0 -s c 0;
    0 0 0 0 0 1;]
k_ele=T\k_local*T;



function k_t=assemPlaneFrame(k_t,k_ele,mode1,mode2)
%assemPlaneFrame This function assembles the element stiffness
%matrix k of the plane frame element with nodes
%i and j into the global stiffness matrix K.
%This function returnes the global stiffness
%matrix K after the element stiffness matrix
%k is assembled.


d(1)=3*node1-2;
d(2)=3*node1-1;
d(3)=3*node1;
d(4)=3*node2-2;
d(5)=3*node2-1;
d(6)=3*node2;
for ii=1:6
    for jj=1:6
        k_t(d(ii),d(jj)=k_t(d(ii),d(jj))+k_ele(ii,jj));
    end
end

%~~~main~~~
clc
clear;
node=[1 -60 0 0;
      2 -60 120 0;
      3 60 120 0;
      4 60 0 0];
ele=[1 1 2;
     2 2 3;
     3 3 4];
E=3e7;
I=[200 100 200];
A=10;
n_ele=length(ele(:,1));
l(1:n_ele)=0;
theta(1:n_ele)=0;

for i=1:n_ele
    l(i)=sqrt((node(ele(i,3),2)-node(ele(i,2),2))^2+(node(ele(i,3),3)-node(ele(i,2),3))^2+(node(ele(i,3),4)-mode(ele(i,2),4))^2);
    theta(i)=atand(node(ele(i,3),3)-node(ele(i,2),3))/(node(ele(i,3),2)-node(ele(i,2),2));
end

dof=length(node(:,1))*3;
f=ones(dof,1)*1e8;
f_loc=zeros(6,1);
u=ones(dof,1)*1e6;
K=zeros(dof);
stress=zeros(n_ele,1);

for i=1:n_ele
    k_ele=PlaneFrameElementStiffness(A,E,I(i),theta(i));
    K=assemPlaneFrame(K,k_ele,ele(i,2),ele(i,3));
end

f(1)=0;
f(2)=0;
f(3)=0;
f(4)=10000;
f(5)=0;
f(6)=0;
f(7)=0;
f(8)=0;
f(9)=5000;
%f(6)=0;
%f(7)=-10000;
%f(8)=0;

u(1)=0;
U(2)=0;
u(3)=0;
u(10)=0;
u(11)=0;
u(12)=0;

index=[];
p=[];
for i=1:dof
    if u(i)~=0
        index=[index,i];
        p=[p;f(i)];
    end
end
u(index)=K(index,index)\p;
f=K*u;

for i=1:n_ele
    u1=[u(3*ele(i,2)-2);u(3*ele(i,2)-1);u(3*ele(i,2));u(3*ele(i,3)-2);u(3*ele(i,3)-2);u(3*ele(i,3)-1);u(3*ele(i,3))];
    f_loc=PlaneFrameElementForece(A,E,I(i),theta(i),u1);
end

function [J,grad,H]=costf_mhe_cstr(d_var,x_rec,u,x_last,sen_noise)
global dt N A C B B0 D D0

Q=1;
R=10^(-5);
pi0=10^(-6);

w_es=d_var(2:N+1);
v_es=zeros(N,1);

x0_es=[d_var(1);x_rec(2,1)+sen_noise(1)];

car=5;%2nd order Carleman
Ai=zeros(car,car,N);
E=zeros(car,car,N);
F=zeros(car,1,N);
GuU=zeros(car,car,N);
I=eye(car);
cx=zeros(car,N+1);
cx0=x0_es;
dummy=[cx0;kron(cx0,cx0)];
cx(:,1)=dummy([1 2 3 5 6]);
for i=1:1:N
    Ai(:,:,i)=A+D*w_es(i)+B*u(i);
    F(:,:,i)=B0*u(i)+D0*w_es(i)+C;
    E(:,:,i)=expm(Ai(:,:,i)*dt);
    GuU(:,:,i)=Ai(:,:,i)\(E(:,:,i)-I);
    cx(:,i+1)=E(:,:,i)*cx(:,i)+GuU(:,:,i)*F(:,:,i);
    x_r=cx(:,i+1);
    dummy=[x_r([1 2]);kron(x_r([1 2]),x_r([1 2]))];
    cx(:,i+1)=dummy([1 2 3 5 6]);
    v_es(i)=x_rec(2,i+1)+sen_noise(i+1)-cx(2,i+1);
end

J=v_es'*Q*v_es+w_es*R*w_es'+(x0_es-x_last)'*pi0*(x0_es-x_last);

% Calculate Gradient=======================================================
CC=zeros(1,car);
CC(2)=1;
dEdw=zeros(car,car,N);dxdw=zeros(car,1,N);

for i=1:N
mE1=dt*D;
mE2=dt^2/2*(D*Ai(:,:,i)+Ai(:,:,i)*D);
mE3=dt^3/6*(D*Ai(:,:,i)^2+Ai(:,:,i)*D*Ai(:,:,i)+Ai(:,:,i)^2*D);
mE4=dt^4/24*(D*Ai(:,:,i)^3+Ai(:,:,i)*D*Ai(:,:,i)^2+Ai(:,:,i)^2*D*Ai(:,:,i)+Ai(:,:,i)^3*D);
% mE5=dt^5/120*(D*Ai(:,:,i)^4+Ai(:,:,i)*D*Ai(:,:,i)^3+Ai(:,:,i)^2*D*Ai(:,:,i)^2+Ai(:,:,i)^3*D*Ai(:,:,i)+Ai(:,:,i)^4*D); 
mE5=0;
dEdw(:,:,i)=mE1+mE2+mE3+mE4+mE5;
end
%dxdx
dxdx=[1,0,2*x0_es(1),x0_es(2),0]';
%dvdx0
dv1dx=-CC*E(:,:,1)*dxdx;
dv2dx=-CC*E(:,:,2)*E(:,:,1)*dxdx;
dv3dx=-CC*E(:,:,3)*E(:,:,2)*E(:,:,1)*dxdx;
%dxdw
for i=1:1:N
    dxdw(:,:,i)=dEdw(:,:,i)*cx(:,i)+Ai(:,:,i)\dEdw(:,:,i)*F(:,:,i)+GuU(:,:,i)*D0-Ai(:,:,i)\D*GuU(:,:,i)*F(:,:,i);
end
%dvdw
dv1dw1=-CC*dxdw(:,:,1);
dv1dw2=0;
dv1dw3=0;
dv2dw1=-CC*E(:,:,2)*dxdw(:,:,1);
dv2dw2=-CC*dxdw(:,:,2);
dv2dw3=0;
dv3dw1=-CC*E(:,:,3)*E(:,:,2)*dxdw(:,:,1);
dv3dw2=-CC*E(:,:,3)*dxdw(:,:,2);
dv3dw3=-CC*dxdw(:,:,3);
% Gradient
dJdx0_1=2*(x_last-x0_es)'*pi0*[1;0]+2*v_es(1)*Q*dv1dx+2*v_es(2)*Q*dv2dx+2*v_es(3)*Q*dv3dx;
dJdw1=2*w_es(1)*R+2*v_es(1)*Q*dv1dw1+2*v_es(2)*Q*dv2dw1+2*v_es(3)*Q*dv3dw1;
dJdw2=2*w_es(2)*R+2*v_es(1)*Q*dv1dw2+2*v_es(2)*Q*dv2dw2+2*v_es(3)*Q*dv3dw2;
dJdw3=2*w_es(3)*R+2*v_es(1)*Q*dv1dw3+2*v_es(2)*Q*dv2dw3+2*v_es(3)*Q*dv3dw3;
grad=[dJdx0_1,dJdw1,dJdw2,dJdw3];

% Calculate Hessian=======================================================
% Elements of each element
%d2xd2x
d2xd2x=[0,0,2,0,0]';
%d2vd2x0
d2v1d2x=-CC*E(:,:,1)*d2xd2x;
d2v2d2x=-CC*E(:,:,2)*E(:,:,1)*d2xd2x;
d2v3d2x=-CC*E(:,:,3)*E(:,:,2)*E(:,:,1)*d2xd2x;
%d2vdxdw
d2v1dxdw1=-CC*dEdw(:,:,1)*dxdx;
d2v2dxdw1=-CC*E(:,:,2)*dEdw(:,:,1)*dxdx;
d2v2dxdw2=-CC*dEdw(:,:,2)*dxdx;
d2v3dxdw1=-CC*E(:,:,3)*E(:,:,2)*dEdw(:,:,1)*dxdx;
d2v3dxdw2=-CC*E(:,:,3)*dEdw(:,:,2)*dxdx;
d2v3dxdw3=-CC*dEdw(:,:,3)*dxdx;

% Elements of Hessian matrix
d2Jd2x=2*pi0+2*Q*(dv1dx^2+dv2dx^2+dv3dx^2)+2*v_es(1)*Q*d2v1d2x+2*v_es(2)*Q*d2v2d2x+2*v_es(3)*Q*d2v3d2x;

d2Jdxdw1=2*Q*(dv1dx*dv1dw1+dv2dx*dv2dw1+dv3dx*dv3dw1)+2*Q*(v_es(1)*d2v1dxdw1+v_es(2)*d2v2dxdw1+v_es(3)*d2v3dxdw1);
d2Jdxdw2=2*Q*(dv1dx*dv1dw2+dv2dx*dv2dw2+dv3dx*dv3dw2)+2*Q*(v_es(2)*d2v2dxdw2+v_es(3)*d2v3dxdw2);
d2Jdxdw3=2*Q*(dv1dx*dv1dw3+dv2dx*dv2dw3+dv3dx*dv3dw3)+2*Q*(v_es(3)*d2v3dxdw3);

% Hessian matrix
H=zeros(4,4);
H(1,1)=d2Jd2x;
H(1,2)=d2Jdxdw1;
H(1,3)=d2Jdxdw2;
H(1,4)=d2Jdxdw3;
H(2,1)=H(1,2);
H(3,1)=H(1,3);
H(4,1)=H(1,4);
end
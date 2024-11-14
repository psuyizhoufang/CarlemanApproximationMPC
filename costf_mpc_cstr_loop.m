% Cost function calculated for CMPC
% Assuming no process noises
function [J,dJdu,H]=costf_mpc_cstr_loop(u,x0)
% load carleman_matrices.mat
global dt A C B B0
dt=0.05;
% Prelocating for speed
Ai=zeros(5,5,8);
F=zeros(5,1,8);
E=zeros(5,5,8);
GuU=zeros(5,5,8);
AGuU=zeros(5,5,8);
xxx=zeros(5,1,8);
Ixxx=zeros(5,1,8);
%==========================================================================
% x_kron
I=eye(5);
for i=1:1:4
    Ai(:,:,i)=A+B*u(i);
    F(:,:,i)=B0*u(i)+C;
    E(:,:,i)=expm(Ai(:,:,i)*dt);
    GuU(:,:,i)=Ai(:,:,i)\(E(:,:,i)-I);
    AGuU(:,:,i)=Ai(:,:,i)\(GuU(:,:,i)-dt*I);
end
for i=5:1:8
    Ai(:,:,i)=Ai(:,:,4);
    F(:,:,i)=F(:,:,4);
    E(:,:,i)=E(:,:,4);
    GuU(:,:,i)=GuU(:,:,4);
    AGuU(:,:,i)=AGuU(:,:,4);
end

x_r=x0;
dummy=[x_r([1 2]);kron(x_r([1 2]),x_r([1 2]))];
xxx0=dummy([1 2 3 5 6]);

xxx(:,:,1)=E(:,:,1)*xxx0+GuU(:,:,1)*F(:,:,1);
x_r=xxx(:,:,1);
dummy=[x_r([1 2]);kron(x_r([1 2]),x_r([1 2]))];
xxx(:,:,1)=dummy([1 2 3 5 6]);

Ixxx(:,:,1)=GuU(:,:,1)*xxx0+AGuU(:,:,1)*F(:,:,1);
Ixxxdt=Ixxx(:,:,1);
for i=2:1:8
    xxx(:,:,i)=E(:,:,i)*xxx(:,:,i-1)+GuU(:,:,i)*F(:,:,i);
    x_r=xxx(:,:,i);
    dummy=[x_r([1 2]);kron(x_r([1 2]),x_r([1 2]))];
    xxx(:,:,i)=dummy([1 2 3 5 6]);
    
    Ixxx(:,:,i)=GuU(:,:,i)*xxx(:,:,i-1)+AGuU(:,:,i)*F(:,:,i);
    Ixxxdt=Ixxxdt+Ixxx(:,:,i);
end
%==========================================================================
% u_kron
    ux=zeros(2,1);
    dux=zeros(2,4);
    for i=1:4
    ux=ux+[u(i);u(i)^2];
    dux(:,i)=[1;2*u(i)];
    end
    ux=ux+4*[u(4);u(4)^2];
%==========================================================================    
% weight
    Q=1;
    R=10^(-5);
   
% Jacobian Matrix
    JA=Q*[0,0,0,0,1];
    JB=R*[0,1];
    J0=R*0;
% Cost Function
   J=J0*8*dt+JA*Ixxxdt+JB*ux*dt;   
% Calculate Gradient=======================================================  
% Sensitivity
dE=zeros(5,5,8);
idE=zeros(5,5,8);
for i=1:4
mE1=dt*B;
mE2=dt^2/2*(B*Ai(:,:,i)+Ai(:,:,i)*B);
mE3=dt^3/6*(B*Ai(:,:,i)^2+Ai(:,:,i)*B*Ai(:,:,i)+Ai(:,:,i)^2*B);
mE4=dt^4/24*(B*Ai(:,:,i)^3+Ai(:,:,i)*B*Ai(:,:,i)^2+Ai(:,:,i)^2*B*Ai(:,:,i)+Ai(:,:,i)^3*B);
dE(:,:,i)=mE1+mE2+mE3+mE4;
idE(:,:,i)=mE1*dt/2+mE2*dt/3+mE3*dt/4+mE4*dt/5;
end

for i=5:1:8
    dE(:,:,i)=dE(:,:,4);
    idE(:,:,i)=idE(:,:,4);
end

dxdx=zeros(5,5,8);
idxdx=zeros(5,5,8);
dxdu=zeros(5,1,8);
idxdu=zeros(5,1,8);

for i=1:8
    dxdx(:,:,i)=E(:,:,i);
    idxdx(:,:,i)=GuU(:,:,i);
end
% i=1;
dxdu(:,:,1)=dE(:,:,1)*xxx0+Ai(:,:,1)\dE(:,:,1)*F(:,:,1)+GuU(:,:,1)*B0-Ai(:,:,1)\B*GuU(:,:,1)*F(:,:,1);
idxdu(:,:,1)=idE(:,:,1)*xxx0+Ai(:,:,1)\idE(:,:,1)*F(:,:,1)+AGuU(:,:,1)*B0-Ai(:,:,1)\B*AGuU(:,:,1)*F(:,:,1);

for i=2:1:8
    dxdu(:,:,i)=dE(:,:,i)*xxx(:,:,i-1)+Ai(:,:,i)\dE(:,:,i)*F(:,:,i)+GuU(:,:,i)*B0-Ai(:,:,i)\B*GuU(:,:,i)*F(:,:,i);
    idxdu(:,:,i)=idE(:,:,i)*xxx(:,:,i-1)+Ai(:,:,i)\idE(:,:,i)*F(:,:,i)+AGuU(:,:,i)*B0-Ai(:,:,i)\B*AGuU(:,:,i)*F(:,:,i);
end

G=zeros(5,1,8);
for kappa = 1:8
  for i = kappa:8
    if i > kappa
      dummy = I;
      for j = kappa+1:i-1
        dummy = dummy*dxdx(:,:,j);
      end
      G(:,:,kappa) = G(:,:,kappa)+idxdx(:,:,i)*dummy*dxdu(:,:,kappa);  
    else
      G(:,:,kappa) = G(:,:,kappa)+idxdu(:,:,i);      
    end    
  end
end
% Special case for u(4)
for i=5:1:8
    G(:,:,4)=G(:,:,4)+G(:,:,i);
end
    dux(:,4)=dux(:,4)*5;

dJdu=zeros(4,1);
for i=1:4
    dJdu(i)=JA*G(:,:,i)+JB*dux(:,i)*dt;
end
% Calculate Hessian =======================================================  
H=zeros(4);
dG=zeros(5,3,8);%dG(:,i,kappa)
for i=1:3
    for kappa=i+1:8
        for k=kappa+1:8
            dummy = I;
            for l = i+1:k-1
                if l==kappa
                    dummy = dummy*dE(:,:,l);
                else
                     dummy = dummy*dxdx(:,:,l);
                end
            end
            dG(:,i,kappa) = dG(:,i,kappa)+idxdx(:,:,k)*dummy*dxdu(:,:,i);  
        end
            dummy = I;
            for j = i+1:kappa-1
                        dummy = dummy*dxdx(:,:,j);
            end
     
        dG(:,i,kappa) = dG(:,i,kappa)+idE(:,:,kappa)*dummy*dxdu(:,:,i);       
           
    end
end

    H(1,1)=0.019;
    H(1,2)=JA*dG(:,1,2);
    H(1,3)=JA*dG(:,1,3);
    H(1,4)=JA*(dG(:,1,4)+dG(:,1,5)+dG(:,1,6)+dG(:,1,7)+dG(:,1,8));
    
    H(2,1)=H(1,2);
    H(2,2)=0.017;
    H(2,3)=JA*dG(:,2,3);
    H(2,4)=JA*(dG(:,2,4)+dG(:,2,5)+dG(:,2,6)+dG(:,2,7)+dG(:,2,8));
    
    
    H(3,1)=H(1,3);
    H(3,2)=H(2,3);
    H(3,3)=0.0147;
    H(3,4)=JA*(dG(:,3,4)+dG(:,3,5)+dG(:,3,6)+dG(:,3,7)+dG(:,3,8));

    H(4,1)=H(1,4);
    H(4,2)=H(2,4);
    H(4,3)=H(3,4);
    H(4,4)=0.1006;
end
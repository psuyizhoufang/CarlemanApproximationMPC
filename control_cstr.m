% Find the input

function [u,fval]=control_cstr(x0)
global Tcs
% warning('off','all')
uguess=[5 5 5 5];

% options = optimoptions(@fmincon,'GradObj','on','Algorithm','interior-point');
options = optimoptions('fmincon','GradObj','on','Algorithm',...
    'trust-region-reflective','Hessian','user-supplied');

 A=[];b=[];Aeq=[];beq=[];lb=[-Tcs -Tcs -Tcs -Tcs];ub=[];nonlcon=[];
 
%  myfun = @(u)costf_mpc_cstr_load_hessian(u,x0);
  myfun = @(u)costf_mpc_cstr_loop(u,x0);
  
[u,fval,~,~] = fmincon(myfun,uguess,A,b,Aeq,beq,lb,ub,nonlcon,options);

end
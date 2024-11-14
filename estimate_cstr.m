% Find optimized input sequence

function [x,fval]=estimate_cstr(x_record,m_input,x_last,sen_noise)
global CAs
% warning('off','all')
uguess=[0,0,0,0];

options = optimoptions('fmincon','GradObj','on','Algorithm',...
    'trust-region-reflective','Hessian','user-supplied');
%,'StepTolerance',1e-12
%,'Hessian','user-supplied'
 A=[];b=[];
 Aeq=[];beq=[];
 lb=[0.08-CAs -10 -10 -10];ub=[0.11-CAs 10 10 10];
 
%  nonlcon=@(d_var)sensor_noise(d_var,x_record);
nonlcon=[];
 
 myfun = @(d_var)costf_mhe_cstr(d_var,x_record,m_input,x_last,sen_noise);
 
[x,fval,~,~] = fmincon(myfun,uguess,A,b,Aeq,beq,lb,ub,nonlcon,options);
end
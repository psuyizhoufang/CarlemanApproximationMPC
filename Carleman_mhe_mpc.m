% Carleman MHE+MPC
% Unstable CSTR example
close all
format compact
clear
% profile on
global q CAf Tf V UA k0 ER deltaH rho Cp Tcs CAs TRs 
global dt N xs us
global A C B B0 D D0

% parameters
q = 100;         % liter/min
CAf = 1;         % mol/liter
Tf = 350;        % K
V = 100;         % liter
UA = 5*10^4;     % J/(min K)
k0 = 7.2*10^10;  % min^-1
ER = 8750;       % K
deltaH = -5*10^4; % J/mol
rho = 1000;      % g/liter
Cp = 0.239;      % J/(g K)
Tcs = 311.1;      % K
CAs = 0.093266127453016; % mol/liter
TRs = 3.850293553820693e+02;% K

xs=[CAs;TRs];
us=Tcs;

dt=0.05;
N=3;
sim_window=2.5;
tp=sim_window/dt;

load real_noise2.mat
load carleman_matrices.mat
% load num_H.mat
w=noise;

w_es=zeros(N,tp);% keep the 1st in the sequence
x0_es=zeros(2,tp);
cx_current=zeros(2,tp+N);
x_record=zeros(2,tp+N+1);
m_input=zeros(1,tp+N);

% load the 1st window
x_record(:,1:1+N)=x_1st_window-xs;
m_input(1:N)=0*ones(1,N);
% Operating time
T=0:dt:sim_window+N*dt;
% Initial condition
x_last=[0;0];
x_real_0=x_1st_window(:,1+N);
% Prelocate Carleman matrices
Ai=zeros(5,5,5);
I=eye(5);
cx=zeros(5,N+1);
cx0=x_last;
dummy=[cx0;kron(cx0,cx0)];
cx(:,1)=dummy([1 2 3 5 6]);

for j=1:tp
% MHE
[design_var,~]=estimate_cstr(x_record(:,j:j+N),m_input(j:j+N-1),x_last,sen_noise(j:j+N));
x0_es(1,j)=design_var(1);
x0_es(2,j)=x_record(2,j)+sen_noise(j);
x_last=x0_es(:,j);

w_es(:,j)=design_var(2:4);
cx0=x0_es(:,j);
dummy=[cx0;kron(cx0,cx0)];
cx(:,1)=dummy([1 2 3 5 6]);
    for k=1:1:N
        Ai(:,:,k)=A+B*m_input(j+k-1)+D*w_es(k,j);
        cx(:,k+1)=expm(Ai(:,:,k)*dt)*cx(:,k)+Ai(:,:,k)\(expm(Ai(:,:,k)*dt)-I)*(C+B0*m_input(j+k-1)+D0*w_es(k,j));
        x_r=cx(:,k+1);
        dummy=[x_r([1 2]);kron(x_r([1 2]),x_r([1 2]))];
        cx(:,k+1)=dummy([1 2 3 5 6]);
    end
    cx_current(:,j+N)=cx(1:2,N+1);

% MPC    
[u,~]=control_cstr(cx_current(:,j+N));
m_input(j+N)=u(1);
[~,x_real]=ode45(@ode_feedback_system,[dt*(j+N),dt*(j+N+1)],x_real_0,odeset,m_input(j+N)+us,noise(j+N));
x_record(:,j+N+1)=x_real(end,:)'-xs;
x_real_0=x_real(end,:);
end

% Simulate Open-loop system for comparison
x0=xs;
x_open=zeros(2,tp+N);
x_open(:,1)=x0;
for j=1:tp+N
    [~,x]=ode45(@ode_open_loop_system,[dt*(j-1),dt*j],x0,odeset,noise(j));
    x_open(:,j+1)=x(end,:);
    x0=x(end,:);
end

figure(1)
subplot(3,1,1);
hold on
box on
plot(T(1:tp+N),x_record(1,1:tp+N)+xs(1),'b-','LineWidth',3)
plot(T(1+N:tp+N),cx_current(1,1+N:tp+N)+xs(1),'r--','LineWidth',3)
axis([0 sim_window+N*dt 0.08 0.11])
xlabel('Time [min]','FontSize', 14)
ylabel('C_A [mol/liter]','FontSize', 14)
grid on

subplot(3,1,2);
hold on
box on
plot(T(1:tp+N),x_record(2,1:tp+N)+xs(2),'b-','LineWidth',3)
plot(T(1+N:tp+N),cx_current(2,1+N:tp+N)+xs(2),'r--','LineWidth',3)
axis([0 sim_window+N*dt 384.8 385.6])
xlabel('Time [min]','FontSize', 14)
ylabel('T_R [K]','FontSize', 14)
grid on

subplot(3,1,3);
hold on
box on
plot(T(1:tp),noise(1:tp),'b-','LineWidth',3)
plot(T(1:tp),w_es(1:tp),'r--','LineWidth',3)
xlabel('Time [min]','FontSize', 14)
ylabel('w [L/min]','FontSize', 14)
grid on

figure(2)
subplot(3,1,1);
hold on
box on
plot(T(1:tp+N),x_record(1,1:tp+N)+xs(1),'b-','LineWidth',3)
plot(T(1:tp+N),x_open(1,1:tp+N),'r-','LineWidth',3)
axis([0 sim_window+N*dt 0.08 0.12])
xlabel('Time [min]','FontSize', 14)
ylabel('C_A [mol/liter]','FontSize', 14)
grid on

subplot(3,1,2);
hold on
box on
plot(T(1:tp+N),x_record(2,1:tp+N)+xs(2),'b-','LineWidth',3)
plot(T(1:tp+N),x_open(2,1:tp+N),'r-','LineWidth',3)
axis([0 sim_window+N*dt 382 387])
xlabel('Time [min]','FontSize', 14)
ylabel('T_R [K]','FontSize', 14)
grid on

subplot(3,1,3);
hold on
box on
plot(T(1:tp+N),m_input+us,'b-','LineWidth',3)
plot(T(1:tp+N),Tcs*ones(1,tp+N),'r-','LineWidth',3)
axis([0 sim_window+N*dt 303 318])
xlabel('Time [min]','FontSize', 14)
ylabel('T_c [K]','FontSize', 14)
grid on

figure(4)
hold on
box on
plot(T(1:tp+N),x_record(1,1:tp+N)+xs(1),'b-','LineWidth',3)
plot(T(1+N:tp+N),cx_current(1,1+N:tp+N)+xs(1),'r--','LineWidth',3)
axis([0 sim_window+N*dt 0.08 0.11])
yticks([0.08 0.09 0.10 0.11])
set(gca,'fontsize',20)
xlabel('Time [min]','FontSize', 20)
ylabel('C_A [mol/liter]','FontSize', 20)
grid on

figure(5)
hold on
box on
plot(T(1:tp+N),x_record(2,1:tp+N)+xs(2),'b-','LineWidth',3)
plot(T(1+N:tp+N),cx_current(2,1+N:tp+N)+xs(2),'r--','LineWidth',3)
axis([0 sim_window+N*dt 384 386])
xlabel('Time [min]','FontSize', 20)
ylabel('T_R [K]','FontSize', 20)
yticks([384 384.5 385 385.5 386])
set(gca,'fontsize',20)
grid on
% ==========Figures for Papers=============================================
figure(11)

hold on
box on
plot(T(1:tp+N),x_record(1,1:tp+N)+xs(1),'b-','LineWidth',3)
% plot(T(1:tp+N),x_open(1,1:tp+N),'r-','LineWidth',3)
axis([0 sim_window+N*dt 0.08 0.11])
set(gca,'fontsize',24)
xlabel('Time [min]','FontSize', 24)
ylabel('C_A [mol/liter]','FontSize', 24)
yticks([0.08 0.09 0.10 0.11])
grid on

figure(12)
hold on
box on
plot(T(1:tp+N),x_record(2,1:tp+N)+xs(2),'b-','LineWidth',3)
% plot(T(1:tp+N),x_open(2,1:tp+N),'r-','LineWidth',3)
axis([0 sim_window+N*dt 383 387])
set(gca,'fontsize',24)
xlabel('Time [min]','FontSize', 24)
ylabel('T_R [K]','FontSize', 24)
yticks([383 384 385 386 387])
grid on

figure(22)
subplot(2,2,1);
hold on
box on
plot(T(1:tp+N),x_record(1,1:tp+N)+xs(1),'b-','LineWidth',3)
plot(T(1:tp+N),x_open(1,1:tp+N),'r-','LineWidth',3)
axis([0 sim_window+N*dt 0.08 0.12])
set(gca,'fontsize',20)
xlabel('Time [min]','FontSize', 20)
ylabel('C_A [mol/liter]','FontSize', 20)
grid on

subplot(2,2,3);
hold on
box on
plot(T(1:tp+N),x_record(2,1:tp+N)+xs(2),'b-','LineWidth',3)
plot(T(1:tp+N),x_open(2,1:tp+N),'r-','LineWidth',3)
axis([0 sim_window+N*dt 382 387])
set(gca,'fontsize',20)
xlabel('Time [min]','FontSize', 20)
ylabel('T_R [K]','FontSize', 20)
grid on

figure(13)
subplot(2,1,1);
hold on
box on
plot(T(1:tp+N),x_record(1,1:tp+N)+xs(1),'b-','LineWidth',3)
% plot(T(1:tp+N),x_open(1,1:tp+N),'r-','LineWidth',3)
axis([0 sim_window+N*dt 0.08 0.11])
set(gca,'fontsize',20)
xlabel('Time (min)','FontSize', 20)
ylabel('C_A (mol/liter)','FontSize', 20)
grid on

subplot(2,1,2);
hold on
box on
plot(T(1:tp+N),x_record(2,1:tp+N)+xs(2),'b-','LineWidth',3)
% plot(T(1:tp+N),x_open(2,1:tp+N),'r-','LineWidth',3)
axis([0 sim_window+N*dt 383 387])
set(gca,'fontsize',20)
xlabel('Time (min)','FontSize', 20)
ylabel('T_R (K)','FontSize', 20)
grid on

figure(23)
subplot(2,1,1);
hold on
box on
plot(T(1:tp+N),x_record(1,1:tp+N)+xs(1),'b-','LineWidth',3)
plot(T(1:tp+N),x_open(1,1:tp+N),'r-','LineWidth',3)
axis([0 sim_window+N*dt 0.08 0.12])
set(gca,'fontsize',20)
xlabel('Time (min)','FontSize', 20)
ylabel('C_A (mol/liter)','FontSize', 20)
grid on

subplot(2,1,2);
hold on
box on
plot(T(1:tp+N),x_record(2,1:tp+N)+xs(2),'b-','LineWidth',3)
plot(T(1:tp+N),x_open(2,1:tp+N),'r-','LineWidth',3)
axis([0 sim_window+N*dt 382 387])
set(gca,'fontsize',20)
xlabel('Time (min)','FontSize', 20)
ylabel('T_R (K)','FontSize', 20)
grid on

figure(24)
hold on
box on
% plot(T(1:tp+N),x_record(1,1:tp+N)+xs(1),'b-','LineWidth',3)
plot(T(1:tp+N),x_open(1,1:tp+N),'r-','LineWidth',3)
axis([0 sim_window+N*dt 0.08 0.12])
yticks([0.08 0.09 0.10 0.11 0.12])
set(gca,'fontsize',20)
xlabel('Time [min]','FontSize', 20)
ylabel('C_A [mol/liter]','FontSize', 20)
grid on

figure(25)
hold on
box on
% plot(T(1:tp+N),x_record(2,1:tp+N)+xs(2),'b-','LineWidth',3)
plot(T(1:tp+N),x_open(2,1:tp+N),'r-','LineWidth',3)
axis([0 sim_window+N*dt 382 387])
yticks([382 383 384 385 386 387])
set(gca,'fontsize',20)
xlabel('Time [min]','FontSize', 20)
ylabel('T_R [K]','FontSize', 20)
grid on
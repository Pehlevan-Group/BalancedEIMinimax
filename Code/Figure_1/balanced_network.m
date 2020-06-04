% spiking model for firing rate predictions/objective
function [rate,re,ri,W_EE,W_EI,W_IE,W_II,F,NE,N,tau] = balanced_network(tau)

% input
% tau: time constant of the network


% returns
% rate: firing rate of input
% re: firing rate of E neurons
% ri: firing rate of I neurons
% W_EE,W_IE,W_II, W_EI, F: weights of the network

%% 
rng(1,'twister');

addpath('../')
set(0,'DefaultAxesColorOrder',brewermap(NaN,'Paired'))
color = get(gca,'colororder');

%%

NE = 20;
NI = 20;
N = NE + NI;

m = 20; % dimension of input

max_rate = 150;  % maximum firing rate in input

% generate weight matrices that satisfy the conditions 
rmatrix = rand(NE, NE);
W_EE = rmatrix'*rmatrix/NE;
[V, D] = eig(W_EE);
alpha = 1./max(abs(diag(D))+1);
W_EI = rand(NE, NI);
W_II = W_EI'*diag(alpha)*W_EI;
W_IE = W_EI';
W_EE = W_EE-eye(NE)-diag(diag(W_EE));
v_th_e = -1/2*diag(W_EE);
v_th_i = 1/2*diag(W_II);

% simulation time 
time = 100*tau;
dt = 1e-4*tau;
time_all = dt:dt:time;

% generate input spike train
x = binornd(1,dt*max_rate*kron(rand(m, 1), ones(1, round(length(time_all)))));

% feed-forward weight
F = normrnd(0.1,0.5,NE,size(x,1));

% initialize the variables
spike_e = zeros(NE,length(time_all)+1);
spike_i = zeros(NI,length(time_all)+1);
y_e = normrnd(0, 1, NE, length(time_all)+1);
y_i = normrnd(0, 1, NI, length(time_all)+1);
rate_e = zeros(NE, length(time_all)+1);
rate_i = zeros(NI, length(time_all)+1);
input = zeros(size(x,1), length(time_all)+1);

for i0 = 1:length(time_all)
    spike_e(:, i0) =  (y_e(:, i0) > v_th_e);
    spike_i(:, i0) =  (y_i(:, i0) > v_th_i);
    input_ee(:, i0) =  dt*(W_EE*rate_e(:,i0) + F*input(:,i0));
    input_ie(:, i0) = - dt*(W_EI*rate_i(:,i0));
    input_ii(:, i0) =  - dt*(W_II*rate_i(:,i0));
    input_ei(:, i0) =  dt*(W_IE*rate_e(:,i0));
    y_e(:, i0 + 1) = y_e(:,i0) + dt/tau*(-y_e(:,i0) + tau/dt*W_EE*spike_e(:,i0) - tau/dt*W_EI*spike_i(:,i0) + tau/dt*F*x(:,i0));
    y_i(:, i0 + 1) = y_i(:,i0) + dt/tau*(-y_i(:,i0) + tau/dt*W_IE*spike_e(:,i0) - tau/dt*W_II*spike_i(:,i0));
    rate_e(:, i0+1) = rate_e(:,i0) + dt/tau*(-rate_e(:,i0) + tau/dt*spike_e(:,i0));
    rate_i(:, i0+1) = rate_i(:,i0) + dt/tau*(-rate_i(:,i0) + tau/dt*spike_i(:,i0));
    input(:,i0+1) = input(:,i0) + dt/tau*(-input(:,i0) + tau/dt*x(:,i0));
end

% average firing rate
rate = mean(input(:,round(2 + length(time_all)/2):round(length(time_all)+1))')'/tau;
re = mean(rate_e(:,round(2 + length(time_all)/2):round(length(time_all)+1))')'/tau;
ri = mean(rate_i(:,round(2 + length(time_all)/2):round(length(time_all)+1))')'/tau;


% plotting the E and I inputs
time_s = (0:5*round(tau/dt))*dt;
start_time = 30000;
plot_time = 5*round(tau/dt);

plot(time_s/tau,input_ee(11,start_time:start_time+plot_time)/max(input_ee(11,start_time:start_time+plot_time)), ...
    'linewidth',2,'linestyle','-','color',color(5,:));hold on
plot(time_s/tau,-input_ie(11,start_time:start_time+plot_time)/max(input_ee(11,start_time:start_time+plot_time)), ...
    'linewidth',2,'linestyle','-','color',color(2,:));hold on

set(gca,'fontsize',15,'linewidth',2)
axis tight
xlabel('Time (\tau)')
ylabel('Input')
title(['\tau=', num2str(tau)]);

legend('E','I')

end
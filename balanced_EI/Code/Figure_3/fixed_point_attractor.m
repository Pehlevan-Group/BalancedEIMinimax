% spiking model with learned weights 
% input: index - the index of the attractor state to recover
%        noise - noise during initialization
function [] = fixed_point_attractor(index, noise)
tau_E = 0.5;
tau_I = 0.4;
load('./Data/optim_w.mat') % load the learned weights

NE = size(r_e, 1);
NI = size(r_i, 1);
N = NE+NI;

W_EE = reshape(w(1:NE^2), NE, NE);
W_II = reshape(w(1+NE^2:NE^2+NI^2), NI, NI);
W_EI = reshape(w(NE^2+NI^2+1:end-2), NI, NE)';
W_IE = W_EI';
W_II = W_II + diag(ones(1, size(W_II,1)));
W_EE = W_EE - diag(ones(1, size(W_II,1)));
thetaE = w(end-1);
thetaI = w(end);
v_th_e = -1/2*diag(W_EE)+thetaE;
v_th_i = 1/2*diag(W_II)+thetaI;

time = 50;
dt = 0.001;
time_all = dt:dt:time;
sig = noise;
spike_e = zeros(NE, length(time_all));
spike_i = zeros(NI, length(time_all));
y_e = thetaE*ones(NE, length(time_all));
y_i = thetaI*ones(NI, length(time_all));
input_e = zeros(NE, length(time_all));
input_i = zeros(NI, length(time_all));

% initialize with some spike rate (poisson spikes) with noise
r_e = (r_e(:,index) + abs(normrnd(0,sig,NE,1)))/tau_E;
r_i = (r_i(:,index) + abs(normrnd(0,sig,NI,1)))/tau_I;
ini_time = 10;
spike_e(:, 1:round(ini_time/dt)) = binornd(1, dt*r_e*ones(1,round(ini_time/dt)));
spike_i(:, 1:round(ini_time/dt)) = binornd(1, dt*r_i*ones(1,round(ini_time/dt)));

for i0 = 2:round(ini_time/dt)
    input_e(:,i0) = input_e(:,i0-1) - dt/tau_E*(input_e(:,i0-1)) + spike_e(:,i0-1);
    input_i(:,i0) = input_i(:,i0-1) - dt/tau_I*(input_i(:,i0-1)) + spike_i(:,i0-1);
end

start_time = round(ini_time/dt) + 1;

% simulation of the spiking dynamics
for i0 = start_time:length(time_all)
    spike_e(:,i0) = (y_e(:,i0-1) > v_th_e);
    spike_i(:,i0) = (y_i(:,i0-1) > v_th_i);
    y_i(:,i0) = W_IE*input_e(:,i0-1) - W_II*input_i(:,i0-1);
    y_e(:,i0) = W_EE*input_e(:,i0-1) - W_EI*input_i(:,i0-1);
    input_e(:,i0) = input_e(:,i0-1) - dt/tau_E*(input_e(:,i0-1)) + spike_e(:,i0-1);
    input_i(:,i0) = input_i(:,i0-1) - dt/tau_I*(input_i(:,i0-1)) + spike_i(:,i0-1);
end

imagesc(time_all, 1:15, input_e)
colormap jet
hold on
line([10,10],[0,15],'linestyle','-','linewidth',4,'color','k')
xlabel('time')
ylabel('neuron')
set(gca,'fontsize',20,'linewidth',2)
end
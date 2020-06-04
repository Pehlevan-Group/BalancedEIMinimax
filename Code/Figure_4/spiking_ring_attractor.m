% spiking model for ring attractor
% input:
%   theta_1 : inhomogeneous input theta_0 during 0-10s (ring attractor)
%   theta_2 : inhomogeneous input theta_0 during 20-30s (ring attractor)
%   gamma : angular frequency for travelling wave

function [] = spiking_ring_attractor(theta_1,theta_2,gamma)
tau = 0.1;

w0 = 0.5;
w1 = 2.7;
h0 = 1/tau;

if gamma == 0 % ring attractor
    h1 = 5;
    time = 40;
else  % travelling wave
    h1 = 0;
    time=20;
end
noise = 0;

NE = 100;
NI = 100;
N = NE+NI;

% weight design
W_IE = 1/sqrt(NE)*sqrt(3*(w0+w1)/NI)*ones(NI, NE);
W_EI = W_IE';
W_II = 2*eye(NI);
cos_mat = zeros(NE, NE);
sin_mat = zeros(NE, NE);
for i = 1:NE
    cos_mat = cos_mat + w1*diag(cos(2*pi*(i-1)/NE)*ones(NE-i+1,1),i-1);
    sin_mat = sin_mat - w1*gamma*diag(sin(2*pi*(i-1)/NE)*ones(NE-i+1,1),i-1);
end
cos_mat = (cos_mat+cos_mat'-diag(diag(cos_mat)));
sin_mat = (sin_mat-sin_mat');
cos_mat = 1/NE*(cos_mat + sin_mat + w0*ones(NE, NE));
W_EE = 1/3*W_EI*W_IE + cos_mat;
thetaI = zeros(NI, 1);
thetaE = zeros(NE, 1);
W_II = W_II + diag(ones(1, size(W_II,1)));
W_EE = W_EE - diag(ones(1, size(W_EE,1)));
v_th_e = -1/2*diag(W_EE)+thetaE;
v_th_i = 1/2*diag(W_II)+thetaI;


dt = 0.001;
time_all = dt:dt:time;

spike_e = zeros(NE, length(time_all));
spike_i = zeros(NI, length(time_all));
y_e = zeros(NE, length(time_all));
y_i = zeros(NI, length(time_all));
input_ee = zeros(NE, length(time_all));
input_ie = zeros(NE, length(time_all));
input_ii = zeros(NI, length(time_all));
input_ei = zeros(NI, length(time_all));
input = zeros(NE, length(time_all));
rate_e = zeros(NE, length(time_all));
rate_i = zeros(NI, length(time_all));
F = eye(NE);

% input
phase = theta_1;
ini_time = 10;
change_time = 20:30;
allphase = 2*pi*(0:NE-1)/NE;
x = binornd(1, dt*h0*ones(NE,round(time/dt)));
x(:, 1:round(ini_time/dt)) = binornd(1, dt*(h0 + h1*cos(phase-allphase)')*ones(1,round(ini_time/dt)));
phase = theta_2;
x(:, (round(change_time(1)/dt)+1):round(change_time(end)/dt)) = binornd(1, dt*(h0 + h1*cos(phase-allphase)')*ones(1,round(ini_time/dt)));

% simulation of the spiking network
for i0 = 1:length(time_all)
    spike_e(:, i0) =  (y_e(:, i0) > v_th_e);
    spike_i(:, i0) =  (y_i(:, i0) > v_th_i);    
    y_e(:, i0 + 1) = y_e(:,i0) + dt/tau*(-y_e(:,i0) + tau/dt*W_EE*spike_e(:,i0) - tau/dt*W_EI*spike_i(:,i0) + tau/dt*F*x(:,i0)) + sqrt(2*dt)*noise*normrnd(0,1,NE,1);
    y_i(:, i0 + 1) = y_i(:,i0) + dt/tau*(-y_i(:,i0) + tau/dt*W_IE*spike_e(:,i0) - tau/dt*W_II*spike_i(:,i0)) + sqrt(2*dt)*noise*normrnd(0,1,NE,1);
    rate_e(:, i0+1) = rate_e(:,i0) + dt/tau*(-rate_e(:,i0) + tau/dt*spike_e(:,i0));
    rate_i(:, i0+1) = rate_i(:,i0) + dt/tau*(-rate_i(:,i0) + tau/dt*spike_i(:,i0));
    input(:,i0+1) = input(:,i0) + dt/tau*(-input(:,i0) + tau/dt*x(:,i0));
end

if gamma == 0
    imagesc(time_all,1:NE,rate_e/tau,[0,250])
    xlabel('time')
    ylabel('neuron')
    hold on
    line([ini_time,ini_time],[0,NE],'linestyle','-','linewidth',4,'color','k')
    hold on
    line([10,10],[0,NE],'linestyle','-','linewidth',4,'color','k')
    hold on
    line([20,20],[0,NE],'linestyle','-','linewidth',4,'color','k')
    hold on
    line([30,30],[0,NE],'linestyle','-','linewidth',4,'color','k')
    set(gca,'fontsize',15,'linewidth',2)
    colormap(jet)
else
    imagesc(time_all,1:NE,rate_e/tau,[0,250])
    xlabel('time')
    ylabel('neuron')
end
end
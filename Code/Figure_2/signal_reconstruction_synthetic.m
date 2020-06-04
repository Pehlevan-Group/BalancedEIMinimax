% spiking model for signal reconstruction

% input: NE - number of excitatory neurons
%        integrator - 1 for time varying signal, 0 for input with constant
%        firing rate

function [] = signal_reconstruction_synthetic(NE,integrator)

set(0,'DefaultAxesColorOrder',brewermap(NaN,'Paired'))
color = get(gca,'colororder');

rng(1,'twister');

if integrator  % input reconstruction for time varying signal
    tau = 0.5;
    lambda = 1/tau;
    NI = 3;
    beta = 10^(-3)/NE;
    
    % designing the weights
    H = pinv(3*rand(NE, NI));
    G = (-H).*(-H>0);
    F = H.*(H>0);
    W_EI = (F+G)';
    W_II = eye(NI);
    W_IE = F+G;
    W_EE = 2*(F'*G+G'*F)-beta*eye(NE);
    v_th_e = -1/2*diag(W_EE);
    v_th_i = 1/2*diag(W_II);
    inp = H;
    
    
    time = 25;
    dt = 1e-4;
    time_all = dt:dt:time;
    
    % the input, square wave
    s = zeros(NI, length(time_all));
    time_1 = 2:8;
    time_2 = 10:16;
    time_3 = 18:24;
    s(1, round(time_1(1)/dt)+1:round(time_1(end)/dt)) = 70; % square wave input
    s(1, round(time_2(1)/dt)+1:round(time_2(end)/dt)) = 30; % square wave input
    s(1, round(time_3(1)/dt)+1:round(time_3(end)/dt)) = 50; % square wave input
    s(2,:) = 10; % constant input
    s(3,:) = 5;  % constant input

    x = zeros(NI, length(time_all)+1);
    spike_e = zeros(NE,length(time_all)+1);
    spike_i = zeros(NI,length(time_all)+1);
    y_e = zeros(NE,length(time_all)+1);
    y_i = zeros(NI, length(time_all)+1);
    rate_e = zeros(NE, length(time_all)+1);
    rate_i = zeros(NI, length(time_all)+1);
    input = zeros(size(x,1), length(time_all)+1);
  
    for i0 = 1:length(time_all)
        x(:, i0+1) = x(:,i0) + dt*(-lambda*x(:, i0) + s(:,i0)); % compute an integration of the input
        spike_e(:, i0) =  (y_e(:, i0) > v_th_e);
        spike_i(:, i0) =  (y_i(:, i0) > v_th_i);
        input_ee(:, i0) =  dt*(W_EE*rate_e(:,i0) + inp'*input(:,i0));
        input_ie(:, i0) = - dt*(W_EI*rate_i(:,i0));
        input_ii(:, i0) =  - dt*(W_II*rate_i(:,i0));
        input_ei(:, i0) =  dt*(W_IE*rate_e(:,i0));
        y_e(:, i0 + 1) = y_e(:,i0) + dt/tau*(-y_e(:,i0) + tau/dt*W_EE*spike_e(:,i0) - tau/dt*W_EI*spike_i(:,i0) + tau*inp'*s(:,i0) + tau*normrnd(0,0.01,NE,1));
        y_i(:, i0 + 1) = y_i(:,i0) + dt/tau*(-y_i(:,i0) + tau/dt*W_IE*spike_e(:,i0) - tau/dt*W_II*spike_i(:,i0) + tau*normrnd(0,0.01,NI,1));
        rate_e(:, i0+1) = rate_e(:,i0) + dt/tau*(-rate_e(:,i0) + tau/dt*spike_e(:,i0));
        rate_i(:, i0+1) = rate_i(:,i0) + dt/tau*(-rate_i(:,i0) + tau/dt*spike_i(:,i0));
    end

    
    figure(1)
    subplot(2,1,1)
    plot(time_all, s(1,:), 'linewidth', 4);
    set(gca,'fontsize',15,'linewidth',3)
    xlabel('Time(s)')
    ylabel('s')
    ylim([0, max(s(1,:)+0.5)]);
    subplot(2,1,2)
    plot(time_all, x(1,1:end-1), 'linewidth', 4, 'color', color(5,:));
    xlabel('Time(s)')
    ylabel('x')
    set(gca,'fontsize',15,'linewidth',3)
    hold on
    plot(time_all, inp(1,:)*rate_e(:,1:end-1), 'linewidth', 2, 'color',color(6,:));
    set(gca,'fontsize',15,'linewidth',3)
    legend('Target','Prediction')

else    % input reconstruction for spike train with constant firing rate
    NI = 10; % 10-dimensional input
    n = 20; % 20 different inputs
    tau = 1;
    time = 2e3*tau;
    dt = 1e-3*tau;
    time_all = dt:dt:time;
    max_rate = 50;
    m = NI;
    r0 = max_rate*rand(m, n);

    x = binornd(1,kron(dt*r0, ones(1, round(length(time_all)/n))));
    
    % designing the weights
    H = pinv(abs(normrnd(0,1,NE, NI)));
    H = 1/NE*H/std(H(:));
    G = (-H).*(-H>0);
    F = H.*(H>0);
    W_EI = (F+G)';
    W_II = eye(NI);
    W_IE = F+G;
    W_EE = 2*(F'*G+G'*F);
    beta = 6e-3/NE;
    W_EE = W_EE - beta*eye(NE);
    v_th_e = -1/2*diag(W_EE);
    v_th_i = 1/2*diag(W_II);
    inp = H;

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
        input(:,i0+1) = input(:,i0) + dt/tau*(-input(:,i0) + tau/dt*x(:,i0));
        y_e(:, i0 + 1) = y_e(:,i0) + dt/tau*(-y_e(:,i0) + tau/dt*W_EE*spike_e(:,i0) - tau/dt*W_EI*spike_i(:,i0) + tau/dt*inp'*x(:,i0));
        y_i(:, i0 + 1) = y_i(:,i0) + dt/tau*(-y_i(:,i0) + tau/dt*W_IE*spike_e(:,i0) - tau/dt*W_II*spike_i(:,i0));
        rate_e(:, i0+1) = rate_e(:,i0) + dt/tau*(-rate_e(:,i0) + tau/dt*spike_e(:,i0));
        rate_i(:, i0+1) = rate_i(:,i0) + dt/tau*(-rate_i(:,i0) + tau/dt*spike_i(:,i0));
    end

for i0 = 1:n
    rate(:,i0) = mean(input(:,round((i0-1)*length(time_all)/n+1e2):round(i0*length(time_all)/n+1))')';
    re(:,i0) = mean(rate_e(:,round((i0-1)*length(time_all)/n+1e2):round(i0*length(time_all)/n+1))')';
    ri(:,i0) = mean(rate_i(:,round((i0-1)*length(time_all)/n+1e2):round(i0*length(time_all)/n+1))')';
end

    figure(1)
    set(0,'DefaultAxesColorOrder',brewermap(NaN,'Paired'))
    color = get(gca,'colororder');
    a = inp*re/tau;b = rate/tau;
    scatter(b(:),a(:),200,'.k')
    set(gca,'fontsize',15,'linewidth',2)
    box on
    x = 0:max_rate;
    y = x;
    hold on
    plot(x,y,'linewidth',2, 'color', color(5,:))
    xlabel('Input')
    ylabel('Reconstruction')
end
end
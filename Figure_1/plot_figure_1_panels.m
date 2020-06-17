
clear 
close all

addpath('../')
set(0,'DefaultAxesColorOrder',brewermap(NaN,'Paired'))
color = get(gca,'colororder');
%% panel_b&c
% plotting the E-I input and net input for balanced network with \tau=0.08, \tau=0.1 and
% \tau = 0.5
tauall = [0.08,0.1,0.5];

for i = 1:3
    load(['./Data/simulation_for_tau',num2str(tauall(i)),'n1.mat'])
    dt = 1e-4*tau;
    tau = tauall(i);
    
    figure(1)
    subplot(3,1,i)
    time_s = (0:5*round(tau/dt))*dt;
    start_time = 1;
    plot_time = 5*round(tau/dt);
    
    plot(time_s/tau,input_ee(start_time:start_time+plot_time)/max(input_ee(start_time:start_time+plot_time)), ...
         'linewidth',2,'linestyle','-','color',color(5,:));hold on
    plot(time_s/tau,-input_ie(start_time:start_time+plot_time)/max(input_ee(start_time:start_time+plot_time)), ...
         'linewidth',2,'linestyle','-','color',color(2,:));hold on
    legend('excitation','inhibition')
    set(gca,'fontsize',15,'linewidth',2)
    axis tight
    xlabel('Time (\tau)')
    ylabel('Input')
    title(['\tau=', num2str(tau)]);
    
    figure(2)
    subplot(3,1,i)
    plot(time_s/tau,input_ee(start_time:start_time+plot_time)/max(input_ee(start_time:start_time+plot_time)) ...
        +input_ie(start_time:start_time+plot_time)/max(input_ee(start_time:start_time+plot_time)),'linewidth',2,'linestyle','-','color',color(7,:));hold on
    plot(time_s/tau,0.5*dt/(max(input_ee(start_time:start_time+plot_time)))*ones(1,length(time_s/tau)),'linewidth',1,'linestyle','--','color','k')
    set(gca,'fontsize',15,'linewidth',2)
    axis tight
    xlabel('Time (\tau)')
    ylabel('Input')
    title(['\tau=', num2str(tau)]);
end

%% panel_d
% an example of spiking pattern for \tau=0.5

figure(3)

load('./Data/simulation_for_tau0.5n1.mat')

dt = 1e-4*tau;
start_time = 1;
plot_time = 5*round(tau/dt);

spike = [spike_e;spike_i];
plotSpikeRaster(logical(spike(:,start_time:start_time+plot_time)),'plottype','vertline','timeperbin',0.0001);
axis off

%% panel_e
% correlation coefficient of the E and the I inputs

figure(4)
load('./Data/data_corr.mat')

semilogx(tau_all, corr,'.-','linewidth',1.5,'markersize',30,'color',color(2,:))
set(gca,'fontsize',15,'linewidth',2)
xlabel('\tau')
ylabel('\rho')

%% panel_h&i
tau_all = [0.01,0.05,0.08,0.1,0.3,0.5,1,5,10];
load('./Data/data_ISI.mat')

figure(5)
histogram(sig_tot{3}./ave_tot{3}, 10,  'facecolor',color(2,:));
xlabel('CV of ISI');
ylabel('Number of Neurons');
set(gca,'fontsize',15,'linewidth',2);

figure(6)
semilogx(tau_all, cv,'.-','linewidth',1.5,'markersize',30,'color',color(2,:))
set(gca,'fontsize',15,'linewidth',2)
xlabel('$\tau$','interpreter','latex')
ylabel('average CV','interpreter','latex')

%% panel f&g

figure(7)

% tau = 0.01
error = 0.1;
load(['./Data/error_', num2str(error), '.mat'])
subplot(1,2,1)
scatter(ropt_e(:)/tau, re(:),40,color(5,:),'filled');
hold on

ylabel('Spiking simulation')
xlabel('Prediction')
set(gca,'fontsize',15,'linewidth',2)
box on
axis tight
hold on
scatter(ropt_i(:)/tau,ri(:),40,color(2,:),'filled');
hold on

plot(re, re, 'linestyle','-','linewidth',2,'color','k');
legend('E','I')
ylabel('Spiking simulation')
xlabel('Prediction')
set(gca,'fontsize',15,'linewidth',2)
box on
axis tight


% tau = 0.1
error = 1;
load(['./Data/error_', num2str(error), '.mat'])
subplot(1,2,2)
scatter(ropt_e(:)/tau, re(:),40, color(5,:),'filled');
hold on

ylabel('Spiking simulation')
xlabel('Prediction')
set(gca,'fontsize',15,'linewidth',2)
box on
axis tight

scatter(ropt_i(:)/tau,ri(:),40,  color(2,:),'filled');
hold on

plot(re, re, 'linestyle','-','linewidth',2,'color','k');

ylabel('Spiking simulation')
xlabel('Prediction')
set(gca,'fontsize',15,'linewidth',2)
box on
legend('E','I')

axis tight

%% panel_g

% prediction error as a function of \tau

tau_all = [0.01,0.05,0.08,0.1,0.3,0.5,1,5,10];
for i0 = 1:length(tau_all)
    load(strcat('./Data/error_',num2str(tau_all(i0)),'.mat'));
    meanval = mean(re(:));
    error_tau(i0) = sum([(ropt_e(:)/tau-re(:)).^2.]);
end
loglog(tau_all, sqrt(error_tau)/20,'.-','markersize',30,'linewidth',1.5,'color',color(2,:))
set(gca,'fontsize',15,'linewidth',2)
xlabel('\tau')
ylabel('$\frac{1}{N}||\mathbf{d}-\mathbf{r}||_2$','interpreter','latex')
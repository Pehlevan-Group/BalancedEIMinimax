
clear 
close all

addpath('../')
%% panel_a
% ring attractor
figure(1)
load('./Data/ring_attractor.mat');
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

% travelling wave with \gamma=0.08
figure(2)
load('./Data/travelling_wave_0.08.mat')
imagesc(time_all(10001:end)-10,1:NE,rate_e(:,10001:end)/tau)
colormap jet
xlabel('time')
ylabel('neuron')
set(gca,'fontsize',15,'linewidth',2)

% travelling wave with \gamma=0.15
figure(3)
load('./Data/travelling_wave_0.15.mat')
imagesc(time_all(10001:end)-10,1:NE,rate_e(:,10001:end)/tau)
colormap jet
xlabel('time')
ylabel('neuron')
set(gca,'fontsize',15,'linewidth',2)

%% panel_b
% ring attractor state is tightly balanced
dt = 0.001;

figure(4)
load('./Data/ring_attractor.mat');
start_time = 20000;
plot_time = 20*round(tau/dt);
time_s = (0:plot_time)*dt;
plot(time_s/tau,input_ee(start_time:start_time+plot_time)/max(input_ee(start_time:start_time+plot_time)),'linewidth',2,'linestyle','-','color',color(8,:));hold on
plot(time_s/tau,input_ei(start_time:start_time+plot_time)/max(input_ee(start_time:start_time+plot_time)),'linewidth',2,'linestyle','-','color',color(2,:));hold on
set(gca,'fontsize',15,'linewidth',2)
xlabel('time(\tau)')
ylabel('Input')


%% panel_c
% grid attractor
figure(5)
load('./Data/grid_attractor.mat');
imagesc(reshape(input_e/tau,63,63));
axis off
colorbar
title('spiking simulation','fontangle','italic')
set(gca,'fontsize',20,'linewidth',2)

%% panel_d
% grid attractor is tightly balanced
figure(6)
plot((1:2000)*dt/tau,-input_ie(10001:12000)/max(input_ee(10001:12000)),'linewidth',2,'linestyle','-','color',color(2,:));hold on
plot((1:2000)*dt/tau,input_ee(10001:12000)/max(input_ee(10001:12000)),'linewidth',2,'linestyle','-','color',color(8,:));hold on

set(gca,'fontsize',15,'linewidth',2)
axis tight
xlabel('Time (\tau)')
ylabel('Input')
title(['\tau=', num2str(tau)]);
set(gca,'fontsize',15,'linewidth',2)
axis tight
xlabel('Time (\tau)')
ylabel('Input')





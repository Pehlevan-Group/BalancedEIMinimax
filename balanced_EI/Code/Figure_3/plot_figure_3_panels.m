
clear 
close all

addpath('../')
set(0,'DefaultAxesColorOrder',brewermap(NaN,'Paired'))
color = get(gca,'colororder');
%% panel_a&b
% fixed point attractors

colormap jet
dt = 0.001;
time_all = dt:dt:50;
figure(1)
subplot(2,1,1)
load('./Data/fixed_point_attractor_sim_input1sig0.mat');
imagesc(time_all, 1:30, [input_e/tau_E;input_i/tau_I], [0,120])
colormap jet
hold on
line([10,10],[0,30],'linestyle','-','linewidth',4,'color','k')
xlabel('time')
ylabel('neuron')
set(gca,'fontsize',15,'linewidth',2)
subplot(2,1,2)
load('./Data/fixed_point_attractor_sim_input1sig10.mat');
imagesc(time_all, 1:30, [input_e/tau_E;input_i/tau_I],[0,120])
colormap jet
%colormap(brewermap([],'YlGnBu'))
hold on
line([10,10],[0,30],'linestyle','-','linewidth',4,'color','k')
xlabel('time')
ylabel('neuron')
set(gca,'fontsize',15,'linewidth',2)
figure(2)
subplot(2,1,1)
load('./Data/fixed_point_attractor_sim_input2sig0.mat');
imagesc(time_all, 1:30, [input_e/tau_E;input_i/tau_I],[0,120])
colormap jet
%colormap(brewermap([],'YlGnBu'))
hold on
line([10,10],[0,30],'linestyle','-','linewidth',4,'color','k')
xlabel('time')
ylabel('neuron')
set(gca,'fontsize',15,'linewidth',2)
subplot(2,1,2)
load('./Data/fixed_point_attractor_sim_input2sig10.mat');
imagesc(time_all, 1:30, [input_e/tau_E;input_i/tau_I],[0,120])
colormap jet
%colormap(brewermap([],'YlGnBu'))
hold on
line([10,10],[0,30],'linestyle','-','linewidth',4,'color','k')
xlabel('time')
ylabel('neuron')
set(gca,'fontsize',15,'linewidth',2)

%% panel_c

set(0,'DefaultAxesColorOrder',brewermap(NaN,'Paired'))
load('./Data/fixed_point_attractor_sim_input2sig0.mat');
dt = 0.001;
tau = tau_E;
figure(3)
start_time = 30000;
plot_time = 30*round(tau/dt);
time_s = (0:plot_time)*dt;
plot(time_s/tau,input_ee(start_time:start_time+plot_time)/max(input_ee(start_time:start_time+plot_time)),'linewidth',2,'linestyle','-','color',color(8,:));hold on
plot(time_s/tau,-input_ei(start_time:start_time+plot_time)/max(input_ee(start_time:start_time+plot_time)),'linewidth',2,'linestyle','-','color',color(2,:));hold on
set(gca,'fontsize',15,'linewidth',2)
xlabel('time(\tau)')
ylabel('Voltage')
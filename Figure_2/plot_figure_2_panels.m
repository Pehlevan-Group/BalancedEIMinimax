
clear 
close all

addpath('../')
set(0,'DefaultAxesColorOrder',brewermap(NaN,'Paired'))
color = get(gca,'colororder');
%% panel_a
% reconstruction with N_E=60, N_I=10
figure(1)
load('./Data/error_vs_N_60seed_1.mat')

scatter(b(:),a(:),200,'.k')
set(gca,'fontsize',15,'linewidth',2)
box on
x = 0:max_rate;
y = x;
hold on
plot(x,y,'linewidth',2, 'color', color(5,:))
xlabel('Input')
ylabel('Reconstruction')

%% panel_b
% reconstruction of time varying inputs
figure(2)
load('./Data/integrator_Ne30_Ni3.mat')

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
plot(time_all, recons, 'linewidth', 2, 'color',color(6,:));
set(gca,'fontsize',15,'linewidth',3)
legend('Target','Prediction')

%% panel_c
% reconstruction of image patches
figure(3)
load('./Data/A169.mat')
colormap(brewermap([],'YlGnBu'))
load(['./Data/patch_reconstruction_50_',num2str(10),'.mat'])
subplot(2,3,1)
imagesc(reshape(a_all{1},13,13),[min(X(:,50)),max(X(:,50))]);
axis off
subplot(2,3,2)
imagesc(reshape(a_all{3},13,13),[min(X(:,50)),max(X(:,50))]);
axis off
subplot(2,3,3)
imagesc(reshape(a_all{5},13,13),[min(X(:,50)),max(X(:,50))]);
axis off

subplot(2,3,1)
title('Step 1','fontangle','italic','fontweight','normal','fontsize',15)
subplot(2,3,2)
title('Step 100','fontangle','italic','fontweight','normal','fontsize',15)
subplot(2,3,3)
title('Step 1000','fontangle','italic','fontweight','normal','fontsize',15)


load(['./Data/patch_reconstruction_50_',num2str(150),'.mat'])
subplot(2,3,4)
imagesc(reshape(a_all{1},13,13),[min(X(:,50)),max(X(:,50))]);
axis off
subplot(2,3,5)
imagesc(reshape(a_all{3},13,13),[min(X(:,50)),max(X(:,50))]);
axis off
subplot(2,3,6)
imagesc(reshape(a_all{5},13,13),[min(X(:,50)),max(X(:,50))]);
axis off
set(gca,'fontsize',20,'linewidth',2);


figure(4)
colormap(brewermap([],'YlGnBu'))
imagesc(reshape(X(:,50),13,13),[min(X(:,50)),max(X(:,50))]);
axis off
title('Target Image Patch','fontangle','italic','fontweight','normal','fontsize',15)


%% panel_d&e

fig = figure(5);
NI_all = 10:20:150;
load('./Data/A169.mat')
for i0 = 1:length(NI_all)
load(['./Data/patch_reconstruction_50_',num2str(NI_all(i0)),'.mat']);error_all(i0) = error(end);var_all(i0) = variance;
end
color = get(gca,'colororder');
set(fig,'defaultAxesColorOrder',[color(2,:); color(8,:)]);
yyaxis left
plot(NI_all,error_all,'linewidth',2)
ylabel('error','fontangle','italic')
yyaxis right
plot(NI_all,var_all,'linewidth',2)
set(gca,'fontsize',20,'linewidth',2);
ylabel('variance','fontangle','italic')
xlabel('N_I','fontangle','italic')
box off
set(0,'DefaultAxesColorOrder',brewermap(NaN,'Paired'))
figure(6)
load(['./Data/patch_reconstruction_50_',num2str(10),'.mat'])
semilogx(error,'linewidth',2);
hold on
load(['./Data/patch_reconstruction_50_',num2str(70),'.mat'])
semilogx(error,'linewidth',2);
hold on
load(['./Data/patch_reconstruction_50_',num2str(150),'.mat'])
semilogx(error,'linewidth',2);
set(gca,'fontsize',20,'linewidth',2);

legend('N^I=10','N^I=70','N^I=150')
xlabel('timestep','fontangle','italic')
ylabel('error','fontangle','italic')
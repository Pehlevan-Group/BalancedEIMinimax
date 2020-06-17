function ropt = fr_prediction(tau,rate,re,ri,W_EE,W_EI,W_IE,W_II,F,NE,N)
% input
% tau: time constant
% rate: input firing rate
% re: output E neuron firing rate
% ri: output I neuron firing rate
% W_EE,W_EI,W_IE,W_II,F: weights of the network
% NE: number of E neurons
% N: total number of neurons

% return
% ropt: firing rate of [E;I] neurons through optimization(solving KKT conditions)


options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluations', 1e6,'OptimalityTolerance',1e-20,'StepTolerance',1e-6);

x = rate*tau;
r0 = [re;ri]*tau;

flag = zeros(1, size(x,2));
for i0 = 1:size(x,2)
functionopt = @(r) sum(((W_EE*r(1:NE)-W_EI*r(NE+1:end)+F*x(:,i0)).*(r(1:NE))).^2) + sum(((W_IE*r(1:NE)-W_II*r(NE+1:end)).*(r(NE+1:end))).^2);
nonlcon = @(r) silentcon(r, x(:,i0), W_EE, W_EI, W_IE, W_II, F, NE);
while flag(i0) < 1
r_0 = r0(:,i0) + rand(N,1);
[ropt(:,i0),~,flag(i0)] = fmincon(functionopt,r_0,[],[],[],[],zeros(1, length(r_0)),[],nonlcon,options);
end
end


% plotting to show the accuracy of prediction
ropt_e = ropt(1:NE,:);
ropt_i = ropt(NE+1:N,:);
figure
scatter(ropt_e(:)/tau, re(:),400, '.');
hold on

ylabel('Spiking simulation')
xlabel('Prediction')
set(gca,'fontsize',15,'linewidth',3)
box on
axis tight

scatter(ropt_i(:)/tau,ri(:),400, '.');
hold on

plot(re, re, 'linestyle','-','linewidth',2,'color','k');

ylabel('Spiking simulation')
xlabel('Prediction')
set(gca,'fontsize',15,'linewidth',3)
box on
legend('E','I')

% condition for silent neurons/active neurons
function [c,ceq] = silentcon(r, x, W_EE, W_EI, W_IE, W_II, F, NE)
    c = [(W_EE*r(1:NE)-W_EI*r(NE+1:end)+F*x);(W_IE*r(1:NE)-W_II*r(NE+1:end))];
    ceq = [((W_EE*r(1:NE)-W_EI*r(NE+1:end)+F*x)).*r(1:NE);((W_IE*r(1:NE)-W_II*r(NE+1:end))).*r(NE+1:end)];
end
end
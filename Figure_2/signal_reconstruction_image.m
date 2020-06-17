% spiking model for input reconstruction of natural image patches

% input: n - the index of the image patch to reconstruct
%        NI - number of inhibitory neurons

function [] = signal_reconstruction_image(n,NI)
load('./Data/A169.mat'); % dictionary learned from http://www.rctn.org/bruno/sparsenet/
rng(1,'twister');

NE = 400; % the number of basis learned  * 2

N = NE+NI;

tau_E = 1;
tau_I = 1;
tau =1;
time = 1*tau_I;
dt =0.001*tau_I;

time_all = dt:dt:time;

r0 = X(:,n); % the n-th image patch

H = [A,-A]; % double the dictionary size for nonnegativity
L = H'*H;
[V_a,D] = eig(L);

[ind,order] = sort(diag(D),'descend');


% design the weights

D = diag(ind(1:NI));
V = V_a(:,order(1:NI));


V = V*sqrt(D)*sqrt(NI);
D = diag(ones(NI,1))/NI; % modification of V and D in order for stability, will not change VDV';

V_p = V.*(V>0);
V_m = V.*(-V>0);

W_EI = (V_p-V_m)*D;
W_II = D ;
W_IE = W_EI';
W_EE = -2*(V_p*D*V_m' + V_m*D*V_p');
beta =  5e-2/NI;
W_EE = W_EE - beta*eye(NE) ;



v_th_e = -1/2*diag(W_EE);
v_th_i = 1/2*diag(W_II);
inp = H;




spike_e = zeros(NE,length(time_all)+1);
spike_i = zeros(NI,length(time_all)+1);
y_e = v_th_e*ones(1, length(time_all)+1); 

input_e = zeros(NE,length(time_all)+1);




for i0 = 1:length(time_all)

    spike_e(:, i0) =  (y_e(:, i0) > v_th_e);
    input_i = zeros(NI,1000+1);
    y_i = zeros(NI,1000+1);

    input_e(:,i0+1) = input_e(:,i0) - dt/tau_E*(input_e(:,i0)) + delta*spike_e(:,i0);
    
    %  we need the inhibitory neuron to update faster for stability
    for j0 = 1:1e3
    spike_i(:, j0) =  (y_i(:, j0) > v_th_i);
    input_i(:,j0+1) = input_i(:,j0) - dt/tau_I*(input_i(:,j0)) + spike_i(:,j0);
    y_i(:,j0+1) = W_IE*input_e(:,i0+1) - W_II*input_i(:,j0+1) ;
    end
  
    y_e(:,i0+1) = W_EE*input_e(:,i0+1) - W_EI*input_i(:,j0+1) + inp'*r0;


    a = inp*mean(input_e(:,1:i0),2);
    error(i0)= mean((r0(:)-a(:)).^2);
    
end

figure
semilogx(time_all, error,'linewidth',2)
xlabel('time')
ylabel('error')
set(gca,'fontsize',20,'linewidth',2)

figure
subplot(1,2,1)
imagesc(reshape(X(:,n),13,13))
axis off
title('Image patch','fontweight','normal')
subplot(1,2,2)
imagesc(reshape(a, 13,13),[min(X(:,n)),max(X(:,n))]);
axis off
title('Reconstruct','fontweight','normal')
end
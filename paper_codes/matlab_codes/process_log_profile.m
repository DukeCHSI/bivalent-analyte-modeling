clear all; close all; clc;
opts = optimset('MaxIter',1e4,'TolFun',1e-8,'TolX',1e-8);

numT = 1020;

save_dir = strcat('profile_figures/');
load('profile_results/BA1_profile_1020.mat');
load('profile_data/BA1_noisy_data_1020.mat');

% SSE(:,end) = [];
sigma = sqrt(6);
p = 1;
profile_method = 2;
df = 14;

N_sample = 20;
N_order = 3;
ydata = output;
% ydata = reshape(ydata,[1400-120,5]);
% Chosen parameters
IC_params = [ka1, ka2, kd1, kd2];
if numT == 2200
    lower_sample = [1.5e3, 4e-5, 3e-3, 1e-12];
    upper_sample = [3e3, 6e-5, 5e-3, 5*1e-5];
else
    lower_sample = [1.5e3, 4e-5, 3e-3, 1e-12];
    upper_sample = [3e3, 6e-5, 5e-3, 5*1e-4];
end
params1 = IC_params;
% [opt_params,fval] = fmincon(@(params) objective(t,ydata,params,concs,R0s,Rmaxs), IC_params,[],[],[],[],LB,[],[], opts);




N = length(t);
nump = length(IC_params);
% Generate the sampled parameters
sampled_params = zeros(nump, N_sample);
for i=1:nump
    sampled_params(i,:) = 10.^(linspace(log10(lower_sample(i)),log10(upper_sample(i)),N_sample));
end
% sampled_params(4,:) = 10.^(linspace(log10(1e-12),log10(UB(4)),N_sample));
% sampled_params(4,1) = 0;
% sampled_params(:,end+1) = IC_params;
% N_sample = N_sample+1;


if profile_method == 1
    demon = length(t)*length(concs)-df;
    fval = fval/demon;
    SSE = SSE;
    
    for i=1:nump
        delta_val = chi2inv(0.95,df);
        threshold(i) = delta_val*fval + fval;
    end
else
    SSE = SSE/sigma^2;
    SSE = length(t)*length(concs)/2*log(2*pi) + 0.5*SSE;
%     fval = fval/sigma^2;
    for i=1:nump
        delta_val = chi2inv(0.95,df);
        threshold(i) = delta_val/2 + min(SSE,[],'all'); %min(SSE(i,:));%
    end
    
end
    
    
    
figure(1)
semilogx(sampled_params(1,:),SSE(1,:),'LineWidth',2)
hold on
yline(threshold(1),'LineWidth',2)
hold off
xlabel('k_{a1}')
ylabel('Negative log likelihood')
set(gca,'FontSize',14)
box on
legend('Profile','Threshold','FontSize',16)
title('Bivalent Analyte Model-1 with Nominal Dissociation Length: k_{a1} Profile')
% ylim([min(SSE(1,:)), threshold(1)])
save_fig = strcat(save_dir,'ka1_profile');
saveas(gcf,save_fig)


figure(2)
semilogx(sampled_params(2,:),SSE(2,:),'LineWidth',2)
hold on
yline(threshold(2),'LineWidth',2)
hold off
xlabel('k_{a2}')
ylabel('Negative log likelihood')
set(gca,'FontSize',14)
box on
legend('Profile','Threshold','FontSize',16)
title('Bivalent Analyte Model-1 with Nominal Dissociation Length: k_{a2} Profile')
% ylim([min(SSE(2,:)), threshold(2)])
save_fig = strcat(save_dir,'ka2_profile');
saveas(gcf,save_fig)


figure(3)
semilogx(sampled_params(3,:),SSE(3,:),'LineWidth',2)
hold on
yline(threshold(3),'LineWidth',2)
hold off
xlabel('k_{d1}')
ylabel('Negative log likelihood')
set(gca,'FontSize',14)
box on
legend('Profile','Threshold','FontSize',16)
title('Bivalent Analyte Model-1 with Nominal Dissociation Length: k_{d1} Profile')
% ylim([min(SSE(3,:)), threshold(3)])
save_fig = strcat(save_dir,'kd1_profile');
saveas(gcf,save_fig)

figure(4)
semilogx(sampled_params(4,:),SSE(4,:),'LineWidth',2)
hold on
yline(threshold(4),'LineWidth',2)
hold off
xlabel('k_{d2}')
ylabel('Negative log likelihood')
set(gca,'FontSize',14)
box on
legend('Profile','Threshold','FontSize',16)
title('Bivalent Analyte Model-1 with Nominal Dissociation Length: k_{d2} Profile')
% ylim([min(SSE(4,:)), threshold(4)])
% xlim([1e-7, 1e-3])


save_fig = strcat(save_dir,'kd2_profile');
saveas(gcf,save_fig)

save(strcat(save_dir,'profile_results_1020.mat'),'SSE','ka1','ka2','kd1','kd2','threshold','sampled_params')


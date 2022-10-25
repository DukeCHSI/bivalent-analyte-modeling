clear all; close all; clc;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
rng(42);

choose_kd2 = 5;
noise_level = 5;

numT = 2200;
dT_asc = 2;
dT_dis = 2;
T_dis = 420;
T_end = 2200;
%% Initialization
% Chosen parameters
concs = [6.25e-8, 1.25e-7, 2.5e-7, 5e-7, 1e-6];
Rmaxs = [606.2377328	440.7803667	322.6427049	240.6212543	191.5538177];
R0s = [42.86566529	76.18425589	118.5703325	168.18267	222.5108977];
ka1 = 2.25e3; %2.99e3;
ka2 = 7e-5; %7.66e-5;
kd1 = 1.60e-3; %4.76e-2;


choose_kd2 = 5;
kd2 = 1e-5;
save_dir = 'generated_data_noise_5/';
save_file_name = strcat(save_dir,'data_5.mat');


t_asc = 120:dT_asc:T_dis;
t_dis = T_dis:dT_dis:T_end; %423:10:4000; %423:5:1400; %423:30:7200;
t = [t_asc, t_dis(2:end)];

y_all = zeros(length(t),3,length(concs));

for i=1:length(concs)    
    R0 = R0s(i);
    Am = concs(i);
    
    y0 = [Rmaxs(i), 0, 0];
    params = [ka1, ka2, kd1, kd2, Am];
    [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,params), t_asc, y0, opts);
    
    y0 = y_asc(end,:);
    params = [0, 0, kd1, kd2, Am];
    [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,params), t_dis, y0, opts);
    
    y_all(:,:,i) = [y_asc; y_dis(2:end,:)];
end


output = zeros(length(t),5);

for i=1:length(concs)
    output(:,i) = y_all(:,2,i) + y_all(:,3,i);
end

max_output = max(max(output));

figure(1)
hold on
for i=1:length(concs)
    %     output(:,i) = R0s(i) + output(:,i) + 10*rand([length(t),1]);%.*output(:,i);
    output(:,i) = R0s(i) + output(:,i) + normrnd(0,noise_level,[length(t),1]);
%     scatter(t, output(:,i))%,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
    
end
plot(t, output)
line([T_dis, T_dis],[min(min(output))-10, max(max(output))+20],'Color','k','LineStyle','--')
input = t;
legend(num2str(concs(1),'%.2e'),num2str(concs(2),'%.2e'),num2str(concs(3),'%.2e'),num2str(concs(4),'%.2e'),num2str(concs(5),'%.2e'))
[~,~] = title('Simulated Data',...
        {strcat(strcat('k_{a1}=',num2str(ka1,'%.2e')), strcat('-----k_{d1}=',num2str(kd1,'%.2e'))),...
         strcat(strcat('k_{a2}=',num2str(ka2,'%.2e')), strcat('-----k_{d2}=',num2str(kd2,'%.2e'))),...
         strcat(strcat('duration\_t_{asc}=',num2str(t_asc(end) - t_asc(1))),      strcat('-----duration\_t_{dis}=',num2str(t_dis(end) - t_dis(1))))
         },'FontSize',17);
set(gca,'FontSize',14)
box on
savefig(strcat(save_dir,'noisy_data'))

% load('data/dataset2.mat')
% 
% Time = dataset{1};
% RU = dataset{2};
% Time = table2array(Time);
% RU = table2array(RU);
% temp_RU = zeros(size(RU));
% temp_Time = zeros(size(Time));
% for i=1:length(concs)
%     temp = [Time(:,i), RU(:,i)];
%     temp(any(isnan(temp), 2), :) = [];
%     rows_to_remove = temp(:,1) < 120;
%     temp(rows_to_remove,:) = [];
%     
%     temp_Time(1:length(temp(:,1)),i) = temp(:,1);
%     temp_RU(1:length(temp(:,2)),i) = temp(:,2);
% end
% 
% RU = temp_RU;
% Time = temp_Time;
% clear temp_RU temp_Time

% for i=1:length(concs)
%     scatter(Time(:,i),RU(:,i),'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
% end

save(save_file_name,'output','y_all', 'concs', 'Rmaxs', 'R0s','ka1','ka2','kd1','kd2','t','t_asc','t_dis');

function dy = bivalent_rhs(t,y,params)
L = y(1);
X1 = y(2);
X2 = y(3);

Am = params(5);

ka1 = params(1);
ka2 = params(2);
kd1 = params(3);
kd2 = params(4);

% ODE equations
dL = -(2*ka1*Am*L - kd1*X1) - (ka2*X1*L - 2*kd2*X2);
dX1 = (2*ka1*Am*L - kd1*X1) - (ka2*X1*L - 2*kd2*X2);
dX2 = ka2*X1*L - 2*kd2*X2;
dy = [dL; dX1; dX2];
end
clear all; close all; clc;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
rng(42);

choose_kd2 = 5;

%% Initialization
% Chosen parameters
concs = [6.25e-08, 1.25e-07, 2.50e-07, 5.00e-07, 1e-6];
Rmaxs = [378.657252718080,293.123182844610,232.536966828669,172.684445954644,135.585008911978];
R0s = [0.180835108773665,9.36648160410458,22.7122970700652,42.0257240967223,61.1465323468423];

ka1 = 2.99e3; %2.99e3;
ka2 = 3e-5; %7.66e-5;
kd1 = 2e-3; %4.76e-2;

kd2_1 = 1e-4;
kd2_2 = 0;
noise_level = 2.25;


t = 120:2:1020;

t_asc = t(1:151);
t_dis = t(151:end);

y_all = zeros(length(t),3,length(concs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Noiseless %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kd2 = kd2_1;


for i=1:length(concs)
    noise = rand([length(t) 1]);
    
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

% max_output = max(max(output));

figure(1)
hold on
for i=1:length(concs)
%     output(:,i) = R0s(i) + output(:,i) + normrnd(0,0.05*max_output,[length(t),1]);
    output(:,i) = R0s(i) + output(:,i);% + normrnd(0,0.05*max_output,[length(t),1]);
    plot(t, output(:,i),'Color',[0 0.4470 0.7410],'LineWidth',4)
end


kd2 = kd2_2;


for i=1:length(concs)
    noise = rand([length(t) 1]);
    
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

% max_output = max(max(output));
for i=1:length(concs)
%     output(:,i) = R0s(i) + output(:,i) + normrnd(0,0.05*max_output,[length(t),1]);
    output(:,i) = R0s(i) + output(:,i);% + normrnd(0,0.05*max_output,[length(t),1]);
    plot(t, output(:,i),'Color',[0.6350 0.0780 0.1840]		,'LineStyle','-.','LineWidth',4)
end
legend(strcat('k_{d2}=',num2str(kd2_1,'%.2e')),'','','','',strcat('k_{d2}=',num2str(kd2_2,'%.2e')))
[~,~] = title('Comparison of Noiseless Simulations','FontSize',17);%;,...
%         {strcat(strcat('k_{a1}=',num2str(ka1,'%.2e')), strcat('-----k_{d1}=',num2str(kd1,'%.2e'))),...
%          strcat(strcat('k_{a2}=',num2str(ka2,'%.2e')), strcat('-----k_{d2}=',num2str(kd2,'%.2e'))),...
%          strcat(strcat('total\_t_{asc}=',num2str(t_asc(end) - t_asc(1))),      strcat('-----total\_t_{dis}=',num2str(t_dis(end) - t_dis(1))))
%          },
%      'FontSize',17);
set(gca,'FontSize',14)
box on
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Noisy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kd2 = kd2_1;


for i=1:length(concs)
    noise = rand([length(t) 1]);
    
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

% max_output = max(max(output));

figure(2)
hold on
for i=1:length(concs)
%     output(:,i) = R0s(i) + output(:,i) + normrnd(0,0.05*max_output,[length(t),1]);
%     output(:,i) = R0s(i) + output(:,i) + normrnd(0,noise_level*max_output,[length(t),1]);
    output(:,i) = R0s(i) + output(:,i) + normrnd(0,noise_level,[length(t),1]);
end
plot(t, output,'Color',[0.6350 0.0780 0.1840],'LineWidth',3)

% sigma = noise_level*max_output;
% save('noisy_data.mat','t','output','ka1','ka2','kd1','kd2','R0s','Rmaxs','concs','t_asc','t_dis','sigma')


kd2 = kd2_2;


for i=1:length(concs)
    noise = rand([length(t) 1]);
    
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
for i=1:length(concs)
%     output(:,i) = R0s(i) + output(:,i) + normrnd(0,0.05*max_output,[length(t),1]);
%     output(:,i) = R0s(i) + output(:,i) + normrnd(0,noise_level*max_output,[length(t),1]);
    output(:,i) = R0s(i) + output(:,i) + normrnd(0,noise_level,[length(t),1]);
    scatter(t, output(:,i),10,'filled', 'MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410])
end
legend(strcat('k_{d2}=',num2str(kd2_1,'%.2e')),'','','','',strcat('k_{d2}=',num2str(kd2_2,'%.2e')))
[~,~] = title('Comparison of Noisy Simulations for Nominal Dissociation Length','FontSize',17);%;,...
%         {strcat(strcat('k_{a1}=',num2str(ka1,'%.2e')), strcat('-----k_{d1}=',num2str(kd1,'%.2e'))),...
%          strcat(strcat('k_{a2}=',num2str(ka2,'%.2e')), strcat('-----k_{d2}=',num2str(kd2,'%.2e'))),...
%          strcat(strcat('total\_t_{asc}=',num2str(t_asc(end) - t_asc(1))),      strcat('-----total\_t_{dis}=',num2str(t_dis(end) - t_dis(1))))
%          },
%      'FontSize',17);
set(gca,'FontSize',14)
box on
hold off



% input = t;
% 
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
% 
% % for i=1:length(concs)
% %     scatter(Time(:,i),RU(:,i))
% % end


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
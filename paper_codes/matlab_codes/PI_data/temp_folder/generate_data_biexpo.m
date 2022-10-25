clear all; close all; clc;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
rng(42);

choose_kd2 = 5;
noise_level = 2.25;

numT = 1370;
dT_asc = 2;
dT_dis = 2.5;

%% Initialization
% Chosen parameters
concs = [3.13e-08, 6.25e-08, 1.25e-07, 2.50e-07, 5.00e-07];
Rmaxs = [378.657252718080,293.123182844610,232.536966828669,172.684445954644,135.585008911978];
R0s = [0.180835108773665,9.36648160410458,22.7122970700652,42.0257240967223,61.1465323468423];
% [2989.77167975456,5.38789094565668e-05,0.00205143980715123,3.85783598541740e-12]
ka1 = 2.99e3; %1e4; %2.99e3;
ka2 = 5.88e-5; %3e-5; %7.66e-5;
kd1 = 7e-4; %4.76e-2;

if choose_kd2 == 12
    kd2 = 3e-12; %1.75e-3;
    save_file_name = strcat('generated_data_N5_biexpo_',num2str(numT),'T/data_12.mat');
elseif choose_kd2 == 5
    kd2 = 3.8e-5;
    save_file_name = strcat('generated_data_N5_biexpo_',num2str(numT),'T/data_5.mat');
elseif choose_kd2 == 4
    kd2 = 3e-4;
    save_file_name = strcat('generated_data_N5_biexpo_',num2str(numT),'T/data_4.mat');
elseif choose_kd2 == 3
    kd2 = 3e-3;
    save_file_name = strcat('generated_data_N5_biexpo_',num2str(numT),'T/data_3.mat');
elseif choose_kd2 == 2
    kd2 = 3e-2;
    save_file_name = strcat('generated_data_N5_biexpo_',num2str(numT),'T/data_2.mat');
else
    kd2 = 3e-1;
    save_file_name = strcat('generated_data_N5_biexpo_',num2str(numT),'T/data_1.mat');
end

t_asc = 120:dT_asc:420;
t_dis = 420:dT_dis:numT; %423:10:4000; %423:5:1400; %423:30:7200;
t = [t_asc, t_dis(2:end)];

y_all = zeros(length(t),2,length(concs));

for i=1:length(concs)    
    R0 = R0s(i);
    Am = concs(i);
    
    y0 = [Rmaxs(i), 0, 0];
    params = [ka1, ka2, kd1, kd2, Am];
    [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,params), t_asc, y0, opts);
    
    y0 = y_asc(end,:);
%     params = [0, 0, kd1, kd2, Am];
%     [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,params), t_dis, y0, opts);
    
    y_dis(:,1) = y0(2)*exp(-kd1*(t_dis-t_dis(1)));
    y_dis(:,2) = y0(3)*exp(-(kd1+kd2)*(t_dis-t_dis(1)));


    y_all(:,:,i) = [y_asc(:,2:3); y_dis(2:end,:)];
end


output = zeros(length(t),5);

for i=1:length(concs)
    output(:,i) = y_all(:,1,i) + y_all(:,2,i);
end

max_output = max(max(output));

figure(1)
hold on
for i=1:length(concs)
    %     output(:,i) = R0s(i) + output(:,i) + 10*rand([length(t),1]);%.*output(:,i);
    output(:,i) = R0s(i) + output(:,i) + normrnd(0,noise_level,[length(t),1]);
    scatter(t, output(:,i),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
end

input = t;

load('data/dataset2.mat')

Time = dataset{1};
RU = dataset{2};
Time = table2array(Time);
RU = table2array(RU);
temp_RU = zeros(size(RU));
temp_Time = zeros(size(Time));
for i=1:length(concs)
    temp = [Time(:,i), RU(:,i)];
    temp(any(isnan(temp), 2), :) = [];
    rows_to_remove = temp(:,1) < 120;
    temp(rows_to_remove,:) = [];
    
    temp_Time(1:length(temp(:,1)),i) = temp(:,1);
    temp_RU(1:length(temp(:,2)),i) = temp(:,2);
end

RU = temp_RU;
Time = temp_Time;
clear temp_RU temp_Time

for i=1:length(concs)
    scatter(Time(:,i),RU(:,i),'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
end

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
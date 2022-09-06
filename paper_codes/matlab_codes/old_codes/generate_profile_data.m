clear all; close all; clc;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
rng(42);

choose_kd2 = 4;
noise_level = 2.25;

numT = 1020; %1370;
dT_asc = 2.5; %2;
dT_dis = 2.5; %2.5;

%% Initialization
% Chosen parameters
concs = [6.25e-08, 1.25e-07, 2.50e-07, 5.00e-07, 1e-6];
Rmaxs = [378.657252718080,293.123182844610,232.536966828669,172.684445954644,135.585008911978];
R0s = [0.180835108773665,9.36648160410458,22.7122970700652,42.0257240967223,61.1465323468423];
% [2989.77167975456,5.38789094565668e-05,0.00205143980715123,3.85783598541740e-12]
ka1 = 1.69e3; %1e4; %2.99e3;
ka2 = 3.83e-5; %3e-5; %7.66e-5;
kd1 = 2.56e-3; %4.76e-2;

if choose_kd2 == 12
    kd2 = 3e-12; %1.75e-3;
    save_file_name = strcat('log_profile_',num2str(numT),'T_new/data_12.mat');
elseif choose_kd2 == 5
    kd2 = 3e-5;
    save_file_name = strcat('log_profile_',num2str(numT),'T_new/data_5.mat');
elseif choose_kd2 == 4
    kd2 = 3e-4;
    save_file_name = strcat('log_profile_',num2str(numT),'T_new/data_4.mat');
elseif choose_kd2 == 3
    kd2 = 3e-3;
    save_file_name = strcat('log_profile_',num2str(numT),'T_new/data_3.mat');
elseif choose_kd2 == 2
    kd2 = 3e-2;
    save_file_name = strcat('log_profile_',num2str(numT),'T_new/data_2.mat');
else
    kd2 = 3e-1;
    save_file_name = strcat('log_profile_',num2str(numT),'T_new/data_1.mat');
end

t_asc = 120:dT_asc:420;
t_dis = 420:dT_dis:numT; %423:10:4000; %423:5:1400; %423:30:7200;
t = [t_asc, t_dis(2:end)];

y_all = zeros(length(t),3,length(concs));

for i=1:length(concs)    
    R0 = R0s(i);
    Am = concs(i);
    
    y0 = [Rmaxs(i), 0, 0];
    params = [ka1, ka2, kd1, kd2, Am];
    [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,params), t_asc, y0, opts);
    
    y0 = y_asc(end,:);
    params = [0, ka2, kd1, kd2, Am];
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
    scatter(t, output(:,i),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
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
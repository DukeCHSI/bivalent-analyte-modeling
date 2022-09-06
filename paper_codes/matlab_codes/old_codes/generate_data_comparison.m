clear all; close all; clc;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
rng(42);

choose_kd2 = 5;
noise_level = 0;

numT = 1020;
dT_asc = 2;
dT_dis = 2;
T_dis = 420;
T_end = 1020;
%% Initialization
% Chosen parameters
concs = 5*[3.13e-08, 6.25e-08, 1.25e-07, 2.50e-07, 5.00e-07];
Rmaxs = [378.657252718080,293.123182844610,232.536966828669,172.684445954644,135.585008911978];
R0s = [0.180835108773665,9.36648160410458,22.7122970700652,42.0257240967223,61.1465323468423];

% ka1 = 2.99e3; %1e4; %2.99e3;
% ka2 = 5.88e-5; %3e-5; %7.66e-5;
% kd1 = 2.05e-3; %4.76e-2;
ka1 = 1.69e3; %1e4; %2.99e3;
ka2 = 3.83e-5; %3e-5; %7.66e-5;
kd1 = 2.56e-3; %4.76e-2;

t_asc = 120:dT_asc:T_dis;
t_dis = T_dis:dT_dis:T_end;
t = [t_asc, t_dis(2:end)];

%%
kd2 = 3e-4;%1e-5;
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
    output(:,i) = R0s(i) + output(:,i) + normrnd(0,noise_level,[length(t),1]);    
end
plot(t, output,'r','LineWidth',4)

%%
kd2 = 3e-7;
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

for i=1:length(concs)
    output(:,i) = R0s(i) + output(:,i) + normrnd(0,noise_level,[length(t),1]);    
end
plot(t, output,'k--','LineWidth',4)


line([T_dis, T_dis],[min(min(output))-10, max(max(output))+20],'Color','k','LineStyle','--')
input = t;
legend('1.57e-7','3.13e-7','6.25e-7','1.25e-6','2.50e-6')
set(gca,'FontSize',14)
box on





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
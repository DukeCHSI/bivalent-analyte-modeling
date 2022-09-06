clear all; close all; clc;
opts = optimset('Display','iter','MaxIter',1e4,'TolFun',1e-8,'TolX',1e-8);

concs = [3.13e-08, 6.25e-08, 1.25e-07, 2.50e-07, 5.00e-07];
baseline = 120;
association = 300;
dissociation = 900;

load("data/dataset6.mat")

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

RUmaxs = max(max(RU))*ones(1,length(concs));
RU0s = RU(1,:);
% Chosen parameters
IC_params = [1e3, 1e-4, 1e-3, 1e-3 RUmaxs, RU0s];%[1e1, 1e-2, 1e-2, 1e-2, RUmaxs, RU0s];
LB = [0, 0, 0, 0, zeros(size(RUmaxs),'like', RUmaxs), zeros(size(RU0s),'like', RU0s)]; %zeros(size(IC_params),'like', IC_params);
UB = [1e6, 1e-1, 1e-1, 1e-1, 1.2*RUmaxs, 1.2*RU0s];%[1e6, 1e-1, 1e-2, 1e-2];

[params,fval] = fmincon(@(params) objective(Time,RU,params,concs,baseline,association,dissociation), IC_params,[],[],[],[],LB,[],[], opts);
run_bivalent(Time,RU,params,concs,baseline,association,dissociation);

fval

function error = objective(t,Ydata,params,concs,baseline,association,dissociation)
ka1 = params(1);
ka2 = params(2);
kd1 = params(3);
kd2 = params(4);
num_pars = 4;


error = 0;
num_data = 0;
for i=1:length(concs)  
    data = [t(:,i), Ydata(:,i)];
    rows_to_remove = data(:,1) == 0;
    data(rows_to_remove,:) = [];
    rows_asc = data(:,1) < baseline + association;
    data_asc = data(rows_asc,:);
    data(rows_asc,:) = []; 
    
    t_asc = data_asc(:,1);
    t_dis = data(:,1);
    t_dis = t_dis - t_dis(1);
    
    R0 = params(i+5+num_pars);
    Am = concs(i);
    
    y0 = [params(i+num_pars), 0, 0];
    ode_params = [ka1, ka2, kd1, kd2, Am];
    [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_asc, y0);
    
    y0 = y_asc(end,:);
%     ode_params = [0, 0, kd1, kd2, Am];
%     [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_dis, y0);
    
    y_pred = R0 + [y_asc(:,2) + y_asc(:,3); y0(2)*exp(-kd1*t_dis) + y0(3)*exp(-(kd1+kd2)*t_dis)];
    y_data = [data_asc(:,2); data(:,2)];
    error = error + sum((y_pred - y_data).^2);
    
    if length(y_data) ~= length(y_pred)
        disp("Lengths are not equal")
    end
    
    num_data = num_data + length(y_data);
end
error = error;%/num_data;
end

function run_bivalent(t,Ydata,params,concs,baseline,association,dissociation)
ka1 = params(1);
ka2 = params(2);
kd1 = params(3);
kd2 = params(4);
num_pars = 4;
figure;
hold on
for i=1:length(concs)  
    data = [t(:,i), Ydata(:,i)];
    rows_to_remove = data(:,1) == 0;
    data(rows_to_remove,:) = [];
    rows_asc = data(:,1) < baseline + association;
    data_asc = data(rows_asc,:);
    data(rows_asc,:) = []; 
    
    t_asc = data_asc(:,1);
    t_dis = data(:,1);
    t_dis = t_dis - t_dis(1);
    
    R0 = params(i+5+num_pars);
    Am = concs(i);
    
    y0 = [params(i+num_pars), 0, 0];
    ode_params = [ka1, ka2, kd1, kd2, Am];
    [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_asc, y0);
    
    y0 = y_asc(end,:);
%     ode_params = [0, 0, kd1, kd2, Am];
%     [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_dis, y0);
    
    y_pred = R0 + [y_asc(:,2) + y_asc(:,3); y0(2)*exp(-kd1*t_dis) + y0(3)*exp(-(kd1+kd2)*t_dis)];
    y_data = [data_asc(:,2); data(:,2)];
    t_dis = data(:,1);
    t_data = [t_asc; t_dis];
    scatter(t_data, y_data)
    plot(t_data, y_pred, 'LineWidth',4)
end
hold off


end

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
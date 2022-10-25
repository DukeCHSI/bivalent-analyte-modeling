clear all; close all; clc;
opts = optimset('Display','iter','MaxIter',1e4,'TolFun',1e-8,'TolX',1e-8);
load("noisy_data.mat");

concs = concs;
data = output;

t_data = t;
t_asc = t_asc;
t_dis = t_dis;

% Chosen parameters
lower_power = [1, -7, -7, -7];
upper_power = [6, -1, -1, -1];
lower_time = [0, 0, 0, 0, 0];
upper_time = [120, 120, 120, 120, 120];

num_run = 20;

best_fval = 1e64;

for irun = 1:num_run
    k0 = 10.^(rand(1,4).*(upper_power - lower_power) + lower_power);
    tstar0 = rand(1,5).*(upper_time - lower_time) + lower_time;

    IC_params = [k0, Rmaxs, tstar0];
%     IC_params = [ka1, ka2, kd1, kd2, Rmaxs, 30, 50, 60, 70, 75];
    
    LB = zeros(size(IC_params),'like', IC_params);

    t_data_truncated = t_data(61:end);
    t_asc_truncated = t_asc(61:end);
    t_dis_truncated = t_dis;
    data_truncated = data(61:end,:);

    [opt_params,fval] = fmincon(@(params) objective(params,data_truncated,t_data_truncated,t_asc_truncated,t_dis_truncated,concs), IC_params,[],[],[],[],LB,[],[], opts);
    
    if fval < best_fval
        best_fval = fval;
        best_params = opt_params;
        best_inits = IC_params;
    end
end


function err = objective(params,data,t,t_asc,t_dis, concs)
ode_params1 = params(1:4);
ode_params = ode_params1;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

y = zeros(length(t),length(concs));



% Solve the ode
for i=1:length(concs)
    Am = concs(i);
    Rmax = params(4+i);
    ode_params(5) = Am;
    
    t_star = params(4+length(concs)+i);

    if t_star >= 1
        t_0 = t_asc(1) - t_star;
        t_pred_asc = t_0:1:(t_asc(1)-1);
        t_asc_temp = [t_pred_asc, t_asc];
        t_asc_temp = t_asc_temp - t_0;
        y0 = [Rmax, 0, 0];
        [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_asc_temp, y0, opts);
        y_asc(1:length(t_pred_asc),:) = [];
    else
        y0 = [Rmax, 0, 0];
        [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_asc, y0, opts);
    end

    y0 = y_asc(end,:);
    ode_params(5) = 0;
    [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_dis, y0, opts);
    
    y(:,i) = [y_asc(:,2) + y_asc(:,3); y_dis(2:end,2) + y_dis(2:end,3)];
    
    ode_params = ode_params1;
end

err = sum(sum((y-data).^2));
end

function y = run_model(params,data,t,t_asc,t_dis, concs)
ode_params1 = params(1:4);
ode_params = ode_params1;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

y = zeros(length(t),length(concs));



% Solve the ode
for i=1:length(concs)
    Am = concs(i);
    Rmax = params(4+i);
    ode_params(5) = Am;
    
    t_star = params(4+length(concs)+i);

    if t_star >= 1
        t_0 = t_asc(1) - t_star;
        t_pred_asc = t_0:1:(t_asc(1)-1);
        t_asc_temp = [t_pred_asc, t_asc];
        t_asc_temp = t_asc_temp - t_0;
        y0 = [Rmax, 0, 0];
        [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_asc_temp, y0, opts);
        y_asc(1:length(t_pred_asc),:) = [];
    else
        y0 = [Rmax, 0, 0];
        [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_asc, y0, opts);
    end

    y0 = y_asc(end,:);
    ode_params(5) = 0;
    [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_dis, y0, opts);
    
    y(:,i) = [y_asc(:,2) + y_asc(:,3); y_dis(2:end,2) + y_dis(2:end,3)];
    
    ode_params = ode_params1;
end

end

%% bivalent right hand side model
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

function a = insert(a, b, idx)
a = [a(1:length(a) < idx), b, a(1:length(a) >= idx)];
end
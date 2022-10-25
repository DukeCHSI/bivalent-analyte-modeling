clear all; close all; clc;

load('data.mat')

[params, fval] = fminsearch(@(params) objective(params, data, t_asc, t_dis, concs1, concs2, Rmaxs, R0s), IC_params, opts);
yOut = run_bivalent(tdata,Ydata,params,concs,R0s,Rmaxs);

function error = objective(params, data, t_asc, t_dis, concs1, concs2, Rmaxs, R0s)

[RU, ~, ~, ~] = run_monovalent(params, t_asc, t_dis, concs1, concs2, Rmaxs, R0s);

error = sum(sum((data-RU).^2));
end


function [RU, RU1, RU2, RU12] = run_monovalent(params, t_asc, t_dis, concs1, concs2, Rmaxs, R0s)
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
t = [t_asc, t_dis];
RU = zeros(length(t),length(concs1));
RU1 = 0;
RU2 = 0;
RU12 = 0;


num_pars = length(params);
params = [params, 0];
params1 = params;
% Solve the ode
for i=1:length(concs1)
    
    R0 = R0s(i);
    A1m = concs1(i);
    
    % Association
    params(num_pars+1) = A1m;
    
    y0 = [Rmaxs(i), 0]; 
    [~, y_asc] = ode15s(@(t,y) monovalent_rhs(t,y,params), t_asc, y0, opts);
    
    y0 = y_asc(end,:);
    params(num_pars+1) = 0;
    [~, y_dis] = ode15s(@(t,y) monovalent_rhs(t,y,params), [t_asc(end), t_dis], y0, opts);
    
    RU(:,i) = R0 + [y_asc(:,2); y_dis(2:end,2)];
    
    params = params1;
end

end


function dy = monovalent_rhs(t,y,params)
L = y(1);
AL = y(2);

Am = params(3);

ka1 = params(1);
kd1 = params(2);

dL = -ka1*Am*L + kd1*AL;
dAL= ka1*Am*L - kd1*AL;

dy = [dL; dAL];
end

function [RU, RU1, RU2, RU12] = run_parallel(params, t_asc, t_dis, concs1, concs2, Rmaxs, R0s)
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
t = [t_asc, t_dis];
RU1 = zeros(length(t),length(concs1));
RU2 = zeros(length(t),length(concs1));
RU12 = zeros(length(t),length(concs1));
RU = zeros(length(t),length(concs1));

num_pars = length(params);
params = [params, 0, 0];
params1 = params;
% Solve the ode
for i=1:length(concs1)
    
    R0 = R0s(i);
    A1m = concs1(i);
    A2m = concs2(i);
    
    % Association
    params(num_pars+1) = A1m;
    params(num_pars+2) = A2m;
    
    y0 = [Rmaxs(i), 0, 0, 0]; 
    [~, y_asc] = ode15s(@(t,y) parallel_rhs(t,y,params), t_asc, y0, opts);
    
    y0 = y_asc(end,:);
    params(num_pars+1) = 0;
    params(num_pars+2) = 0;
    [~, y_dis] = ode15s(@(t,y) parallel_rhs(t,y,params), [t_asc(end), t_dis], y0, opts);
    
    RU(:,i) = R0 + [y_asc(:,2) + y_asc(:,3) + 2*y_asc(:,4); y_dis(2:end,2) + y_dis(2:end,3) + 2*y_dis(2:end,4)];
    RU1(:,i) = [y_asc(:,2); y_dis(2:end,2)];
    RU2(:,i) = [y_asc(:,3); y_dis(2:end,3)];
    RU12(:,i) = [y_asc(:,4); y_dis(2:end,4)];
    
    params = params1;
end

end

function dy = parallel_rhs(t,y,params)
L = y(1);
A1L = y(2);
A2L = y(3);
A1LA2 = y(4);

A1m = params(9);
A2m = params(10);

ka1 = params(1);
kd1 = params(2);
ka2 = params(3);
kd2 = params(4);
ka3 = params(5);
kd3 = params(6);
ka4 = params(7);
kd4 = params(8);

% ODE equations
dL = - ka1*A1m*L + kd1*A1L - ka2*A2m*L + kd2*A2L;
dA1L = ka1*A1m*L - kd1*A1L - ka4*A1L*A2m + kd4*A1LA2;
dA2L = ka2*A2m*L - kd2*A2L - ka3*A2L*A1m + kd3*A1LA2;
dA1LA2 = ka3*A2L*A1m - kd3*A1LA2 + ka4*A1L*A2m - kd4*A1LA2;

dy = [dL; dA1L; dA2L; dA1LA2];
end
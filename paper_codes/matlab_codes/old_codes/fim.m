clear all; close all; clc;

data = load('generated_data_N5_2200T/data_5.mat');


%% Initialization
t_asc = 120:2:420;
t_dis = 420:2:1400; %423:10:4000; %423:5:1400; %423:30:7200;
t = [t_asc, t_dis(2:end)];
concs = data.concs;
R0s = data.R0s;
Rmaxs = data.Rmaxs;
params = [data.ka1, data.ka2, data.kd1, data.kd2];

% Pertubed percentage
h = 1e-4;

% Intialize the pertubed parameters
params1 = params;
params2 = params;

% Get the number of parameters and the number of data
pars_change = [1, 2, 3, 4];
num_par = length(pars_change);
num_ins = length(t);

% Initialize the two sensitivity matrix
% sum_S: sensitivity matrix for X1 + X2
% sep_S: sensitivity matrix for [X1, X2]
sum_S = zeros(num_ins*length(concs), num_par);

%% Compute the sensitivity matrices

for ipar=1:num_par
    i = pars_change(ipar);
    % Pertube each parameter by 1% above and 1% below
    params1(i) = h*params(i);
    params2(i) = 0.1*params(i);
    
    % Solve the ode with the new parameters
    out = run_bivalent(t,params,concs,R0s,Rmaxs);
    out2 = run_bivalent(t,params2,concs,R0s,Rmaxs);
    
    % Calculate the change of sum using central finite difference
    sum_S(:,ipar) = (out-out2)/(params(i) - params2(i));
    
    % Reset parameters
    params1 = params;
    params2 = params;
end
%% Compute the number of identifiable parameters for X1 + X2
% Remove the outputs at t=1
sum_S = sum_S(2:end,:);
% Compute Fisher Information Matrix
sum_FIM = sum_S'*sum_S;
% Compute the Cramer-Rao bound covariance matrix
sum_C = inv(sum_FIM);
% Check if FIM*C ~ I
invAA = sum_C*sum_FIM
% Compute the rank of FIM
sum_rank = rank(sum_FIM);

% Get standard deviation for parameters
sum_SE = sqrt(diag(sum_C));

% Compute the coefficient of variation
sum_cv = zeros(1,num_par);
for i=1:num_par
    sum_cv(i) = 100*sum_SE(i)/abs(params(i));
end


%%
function y_all = run_bivalent(t,params,concs,R0s,Rmaxs)
t_asc = 120:2:420;
t_dis = 420:2:1400; %423:10:4000; %423:5:1400; %423:30:7200;

ka1 = params(1);
ka2 = params(2);
kd1 = params(3);
kd2 = params(4);

y_all = zeros(length(t),length(concs));

for i=1:length(concs)
    R0 = R0s(i);
    Am = concs(i);
    
    y0 = [Rmaxs(i), 0, 0];
    ode_params = [ka1, ka2, kd1, kd2, Am];
    [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_asc, y0);
    
    y0 = y_asc(end,:);
    ode_params = [0, 0, kd1, kd2, Am];
    [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_dis, y0);
    
    y_all(:,i) = R0 + [y_asc(:,2) + y_asc(:,3); y_dis(2:end,2) + y_dis(2:end,3)];
end

y_all  = y_all(:);
end

%%
function dy = bivalent_rhs(t,y,params)
L  = y(1);
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
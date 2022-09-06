clear all; close all; clc;
opts = optimset('MaxIter',1e4,'TolFun',1e-8,'TolX',1e-8);
numT = 2200;
load(strcat('profile_data/BA1_noisy_data_2200.mat'));
rng(42);
%% Initialization
N_sample = 20;
N_order = 3;

ydata = output;


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

LB = 0*IC_params;
LB1 = LB;

N = length(t);
nump = length(IC_params);
% Generate the sampled parameters
sampled_params = zeros(nump, N_sample);
for i=1:nump
    sampled_params(i,:) = 10.^(linspace(log10(lower_sample(i)),log10(upper_sample(i)),N_sample));
end

% N_sample = N_sample+1;
% sampled_params(:,end+1) = IC_params;
% for i=1:nump
%     sampled_params(i,:) = sort(sampled_params(i,:));
% end

SSE = zeros(nump,N_sample);
for ip=1:nump
    for is=1:N_sample
        formatSpec = 'IP is %f, IS is %f \n';
        fprintf(formatSpec,ip,is)
        
            
            IC_params(ip) = [];
            LB(ip) = [];
            fixed_param = sampled_params(ip,is);
            
    %         try 
                [opt_params,fval] = fmincon(@(params) objective(params,ydata, t_asc, t_dis,ip, fixed_param, concs, Rmaxs, t_stars, R0s, numT), IC_params,[],[],[],[],LB,[],[], opts);

                SSE(ip,is) = fval;
    %         catch
    %             SSE(ip,is) = NaN;
    %         end
            IC_params = params1;
            LB = LB1;
    end
end

save(strcat('profile_results/BA1_profile_2200.mat'),'SSE')

function err = objective(params,data,t_asc, t_dis,ip,fixed_param, concs, Rmaxs, t_stars, R0s, numT)
params1 = insert(params, fixed_param, ip);
t = [t_asc, t_dis(2:end)];


ode_params1 = params1(1:4);
ode_params = ode_params1;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

y = zeros(length(t),length(concs));

% Solve the ode
for i=1:length(concs)
    R0 = R0s(i);
    Am = concs(i);
    Rmax = Rmaxs(i);
    ode_params(5) = Am;
    t_star = t_stars(i);
    
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
    ode_params(1) = 0;
    ode_params(2) = 0;
    ode_params(5) = 0;
    [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_dis, y0, opts);
    
    y(:,i) = [y_asc(:,2) + y_asc(:,3); y_dis(2:end,2) + y_dis(2:end,3)];
    
    ode_params = ode_params1;
end

err = sum(sum((y-data).^2));
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
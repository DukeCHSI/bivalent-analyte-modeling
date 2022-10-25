clear all; close all; clc;
opts = optimset('MaxIter',1e4,'TolFun',1e-8,'TolX',1e-8);
rng(42);


%% Monovalent Fitting
load('generated_data_N5_2200T_monofirst/data_5.mat');
%% Initialization
% Chosen parameters
concs = [3.13e-08, 6.25e-08, 1.25e-07, 2.50e-07, 5.00e-07];
%Rmaxs = [2.53e2, 2.23e2, 1.87e2, 1.56e2, 1.38e2];
%R0s = [5.81, 1.30e1, 2.34e1, 3.83e1, 5.42e1];
Rmaxs = [378.657252718080,293.123182844610,232.536966828669,172.684445954644,135.585008911978];
R0s = [0.180835108773665,9.36648160410458,22.7122970700652,42.0257240967223,61.1465323468423];

true_params = [ka1, ka2, kd1, kd2];

Ydata = output;
tdata = t;

RUmaxs = max(max(Ydata))*ones(1,length(concs));
RU0s = Ydata(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IC_params = [1e3, 1e-3, RUmaxs, RU0s];

% Chosen parameters
LB = zeros(size(IC_params),'like', IC_params);
UB = [1e6, 1e-1, 1e-2 , 1e-2];%[1e6, 1e-1, 1e-2, 1e-2];


[params,fval] = fmincon(@(params) monovalent_objective(t,Ydata,params,concs,R0s,Rmaxs), IC_params,[],[],[],[],LB,[],[], opts);

%% Bivalent Fitting
choose_kd2s = 5;
choose_ICs = 1:6;

% ka1_order = 3:1:5;
ka1_order = log10(params(1));
ka2_order = -5:1:-3;
kd1_order = -5:1:-3;
% kd1_order = log10(params(2));
kd2_order = -6:1:-4;

IC_params_list = [];
for ika1 = 1:length(ka1_order)
    for ika2 = 1:length(ka2_order)
        for ikd1 = 1:length(kd1_order)
            for ikd2 = 1:length(kd2_order)
                IC_params_list = [IC_params_list; ...
                    ka1_order(ika1), ka2_order(ika2), kd1_order(ikd1), kd2_order(ikd2)];
            end
        end
    end
end
IC_params_list = 10.^IC_params_list;


for choose_kd2 = choose_kd2s
    for iIC = 1:length(IC_params_list)
        load('generated_data_N5_2200T_monofirst/data_5.mat');
        kd2_dir = 'generated_data_N5_2200T_monofirst/result_5/';

        %% Initialization
        % Chosen parameters
        concs = [3.13e-08, 6.25e-08, 1.25e-07, 2.50e-07, 5.00e-07];
        %Rmaxs = [2.53e2, 2.23e2, 1.87e2, 1.56e2, 1.38e2];
        %R0s = [5.81, 1.30e1, 2.34e1, 3.83e1, 5.42e1];
        Rmaxs = [378.657252718080,293.123182844610,232.536966828669,172.684445954644,135.585008911978];
        R0s = [0.180835108773665,9.36648160410458,22.7122970700652,42.0257240967223,61.1465323468423];
        
        true_params = [ka1, ka2, kd1, kd2];
        
        Ydata = output;
        tdata = t;
        
        RUmaxs = max(max(Ydata))*ones(1,length(concs));
        RU0s = Ydata(1,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IC_params = [IC_params_list(iIC,:), RUmaxs, RU0s];
        
        IC_name = strcat('IC_',num2str(iIC));
        save_dir = strcat(kd2_dir,IC_name,'/');
        if ~exist(save_dir, 'dir')
            mkdir(save_dir)
        end
        save_fig = strcat(save_dir,'fitting.pdf');
        save_result = strcat(save_dir,'result.mat');
        
        % Chosen parameters
        LB = zeros(size(IC_params),'like', IC_params);
        UB = [1e6, 1e-1, 1e-2 , 1e-2];%[1e6, 1e-1, 1e-2, 1e-2];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        [params,fval] = fmincon(@(params) objective(t,Ydata,params,concs,R0s,Rmaxs), IC_params,[],[],[],[],LB,[],[], opts);
        yOut = run_bivalent(t,Ydata,params,concs,R0s,Rmaxs);
        
        figure(1)
        plot(t,Ydata)
        hold on
        plot(t,yOut,'LineWidth',2)
        hold off
        
        saveas(gcf,save_fig)
        save(save_result,'params','fval','IC_params','true_params','Rmaxs','R0s')
        close all;
    end
end

function error = monovalent_objective(t,Ydata,params,concs,R0s,Rmaxs)
t_asc = 120:2:420;
t_dis = 420:10:2200; %423:10:4000; %423:5:1400; %423:30:7200;
t = [t_asc, t_dis];

ka1 = params(1);
kd1 = params(2);


y_all = zeros(length(t)-1,length(concs));

for i=1:length(concs)
    R0 = params(i+5+2);
    Am = concs(i);
    
    Rmax = params(i+2);
    
    y_asc = Rmax*(1-exp(-(ka1*Am + kd1)*t_asc))/(1+kd1/(ka1*Am));
    
    y0 = y_asc(end);
    y_dis = y0*exp(-kd1*(t_dis-t_dis(1)));
    
    y_all(:,i) = R0 + [y_asc'; y_dis(2:end)'];
end

error = sum(sum((Ydata - y_all).^2));
end



function error = objective(t,Ydata,params,concs,R0s,Rmaxs)
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

t_asc = 120:2:420;
t_dis = 420:10:2200; %423:10:4000; %423:5:1400; %423:30:7200;
t = [t_asc, t_dis];

ka1 = params(1);
ka2 = params(2);
kd1 = params(3);
kd2 = params(4);

y_all = zeros(length(t)-1,length(concs));

for i=1:length(concs)
    R0 = params(i+5+4);
    Am = concs(i);
    
    y0 = [params(i+4), 0, 0];
    ode_params = [ka1, ka2, kd1, kd2, Am];
    [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_asc, y0,options);
    
    y0 = y_asc(end,:);
    ode_params = [0, 0, kd1, kd2, Am];
    [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_dis, y0,options);
    
    y_all(:,i) = R0 + [y_asc(:,2) + y_asc(:,3); y_dis(2:end,2) + y_dis(2:end,3)];
end

error = sum(sum((Ydata - y_all).^2));
end

function y_all = run_bivalent(t,Ydata,params,concs,R0s,Rmaxs)
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

t_asc = 120:2:420;
t_dis = 420:10:2200; %423:10:4000; %423:5:1400; %423:30:7200;
t = [t_asc, t_dis];

ka1 = params(1);
ka2 = params(2);
kd1 = params(3);
kd2 = params(4);

y_all = zeros(length(t)-1,length(concs));

for i=1:length(concs)
    R0 = params(i+5+4);
    Am = concs(i);
    
    y0 = [params(i+4), 0, 0];
    ode_params = [ka1, ka2, kd1, kd2, Am];
    [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_asc, y0,options);
    
    y0 = y_asc(end,:);
    ode_params = [0, 0, kd1, kd2, Am];
    [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_dis, y0,options);
    
    y_all(:,i) = R0 + [y_asc(:,2) + y_asc(:,3); y_dis(2:end,2) + y_dis(2:end,3)];
end
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
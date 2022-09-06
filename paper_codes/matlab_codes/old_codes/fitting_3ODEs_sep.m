clear all; close all;
opts = optimset('Display','iter','MaxIter',10000,'TolFun',1e-12,'TolX',1e-12);
table = readtable('../digitized data/digitized data example 1 newformat v2.xlsx');
clc;

T(:,1) = table.Ligand1_Concentration1;
Y(:,1) = table.Var2;
T(:,2) = table.Ligand1_Concentration2;
Y(:,2) = table.Var4;
T(:,3) = table.Ligand1_Concentration3;
Y(:,3) = table.Var6;
T(:,4) = table.Ligand1_Concentration4;
Y(:,4) = table.Var8;
T(:,5) = table.Ligand1_Concentration5;
Y(:,5) = table.Var10;


Y_temp = Y;
analyte_concs = [6.67e-09 1.33e-08 2.50e-08 5.00e-08 3.33e-07];
temp_analyte_concs = analyte_concs;
RUmax_global = max(max(Y));
RUmax_local = max(Y);


params_save = zeros(length(analyte_concs),5);

for ii=5
    analyte_concs = temp_analyte_concs(ii);
    Y = Y_temp(:,ii);
    
    num_conc = length(analyte_concs);
    
    t = T(:,1);
    
    t_asc = t(1:300);
    Y_asc = Y(1:300,:);
    t_dis = t(301:end);
    Y_dis = Y(301:end,:);
    
    t_dis_start = t(300);
    Y_dis_start = Y(300,:);
    
    dis_start_idx = 301;
    
    RU0 = Y(1,:);
    
    if num_conc == 1
        % Chosen parameters
        IC_params = [1e4, 1e-2, 1e-4, 1e-4, 0.5*RUmax_global];%, 0.1*RUmax, Es, Fs];
        LB = zeros(size(IC_params),'like', IC_params);
        %LB(5:end) = -Inf;
        UB = [1e6, 1e-1, 1e-2, 1e-2, RUmax_global];
    else
        % Chosen parameters
        IC_params = [1e4, 1e-2, 1e-4, 1e-4, 0.5*RUmax_local];%, 0.1*RUmax, Es, Fs];
        LB = zeros(size(IC_params),'like', IC_params);
        %LB(5:end) = -Inf;
        UB = [1e6, 1e-1, 1e-2, 1e-2, RUmax_local];
    end
    [params,fval] = fmincon(@(params) bivalent_obj(t,Y,params,analyte_concs,dis_start_idx,Y_dis_start), IC_params,[],[],[],[],LB,UB,[], opts);
    
    y = run_bivalent(t,Y,params,analyte_concs,dis_start_idx,Y_dis_start);
    
    figure(ii)
    hold on
    for i=1:num_conc
        scatter(t,Y(:,i),'MarkerFaceColor','b')
        plot(t,y(:,i,2)+y(:,i,3),'b','LineWidth',4)
        plot(t,y(:,i,2),'LineWidth',4)
        plot(t,y(:,i,3),'LineWidth',4)
    end
    params_save(ii,:) = params;
end
hold off


save('params_replicate1.mat','params_save')


function y = run_bivalent(t,Ydata,params,analyte_concs,dis_start_idx,Y_dis_start)
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
num_conc = length(analyte_concs);
t_asc = t(1:dis_start_idx-1);
t_dis = t(dis_start_idx:end);
y = zeros(length(t),num_conc,3);

for i=1:num_conc
    if num_conc == 1
        y0 = [params(5), 0, 0];
    else
        y0 = [params(4+i), 0, 0];
    end
    ode_params = params;
    % Association
    [~,yout] = ode15s(@(t,y) bivalent_rhs3(t,y,ode_params,analyte_concs(i)), t_asc, y0, opts);
    y_asc = yout;
    
    % Dissociation
    y0 = yout(end,:);
    ode_params(1) = 0;
    ode_params(2) = 0;
    [~,yout] = ode15s(@(t,y) bivalent_rhs3(t,y,ode_params,analyte_concs(i)), t_dis, y0, opts);
    y_dis = yout;
    
    yout = [y_asc; y_dis];
    y(:,i,:) = yout;
end

end

function err = bivalent_obj(t,Ydata,params,analyte_concs,dis_start_idx,Y_dis_start)
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
num_conc = length(analyte_concs);
err = zeros(length(analyte_concs),1);

t_asc = t(1:dis_start_idx-1);
t_dis = t(dis_start_idx:end);

for i=1:num_conc
    if num_conc == 1
        y0 = [params(5), 0, 0];
    else
        y0 = [params(4+i), 0, 0];
    end
    
    ode_params = params;
    % Association
    [~,yout] = ode15s(@(t,y) bivalent_rhs3(t,y,ode_params,analyte_concs(i)), t_asc, y0, opts);
    y_asc = yout(:,2) + yout(:,3);
    
    % Dissociation
    y0 = yout(end,:);
    ode_params(1) = 0;
    ode_params(2) = 0;
    [~,yout] = ode15s(@(t,y) bivalent_rhs3(t,y,ode_params,analyte_concs(i)), t_dis, y0, opts);
    y_dis = yout(:,2) + yout(:,3);
    
    y = [y_asc; y_dis];
    err(i) = sum(y-Ydata(:,i)).^2;
end

err = sum(err);
end

%%
function dy = bivalent_rhs3(t,y,params,analyte_conc)
L = y(1);
X1 = y(2);
X2 = y(3);

Am = analyte_conc;

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
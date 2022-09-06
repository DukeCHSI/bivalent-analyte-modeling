clear all; close all; clc;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
rng(42);

choose_kd2 = 5;

%% Initialization
% Chosen parameters
concs = [6.25e-08, 1.25e-07, 2.50e-07, 5.00e-07, 1e-6];
Rmaxs = [9.24e+02 7.75e+02 7.02e+02 6.90e+02 7.29e+02];
R0s = 0*[41.7, 73.6, 118, 174, 233];
AL20 = [9.46e+01, 1.04e+02, 1.25e+02, 1.75e+02, 2.64e+02];


ka1 = 1.60e+03; 
ka2 = 6.11e-05;
kd1 = 6.66e-03;
kd2_1 = 1.00e-04;
kd2_2 = 0;

noise_level = 2.25;

% t = 0:2:2200;
% t_asc = t(1:211);
% t_dis = t(211:end);
t = 120:2:2200;
t_asc = t(1:151);
t_dis = t(151:end);

t_stars = [2.75e+02 3.16e+02 3.03e+02 2.46e+02 1.66e+02];

y_all = zeros(length(t),3,length(concs));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Noisy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kd2 = kd2_1;


for i=1:length(concs)
    t_star = t_stars(i);
    
    
    R0 = R0s(i);
    Am = concs(i);
    
    t_0 = t_asc(1) - t_star;
    t_pred_asc = t_0:1:(t_asc(1)-1);
    t_asc_temp = [t_pred_asc, t_asc];
    t_asc_temp = t_asc_temp - t_0;

    y0 = [Rmaxs(i), 0, 0];
    params = [ka1, ka2, kd1, kd2, Am];
    [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,params), t_asc_temp, y0, opts);
    y_asc(1:length(t_pred_asc),:) = [];
    

    y0 = y_asc(end,:);
    params = [ka1, ka2, kd1, kd2, 0];
    [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,params), t_dis, y0, opts);
    
    y_all(:,:,i) = [y_asc; y_dis(2:end,:)];

end





output = zeros(length(t),5);



for i=1:length(concs)
    output(:,i) = y_all(:,2,i) + y_all(:,3,i);
end

figure(2)
hold on
for i=1:length(concs)
    output(:,i) = output(:,i) + normrnd(0,noise_level,[length(t),1]);
end
plot(t, output,'Color',[0.6350 0.0780 0.1840],'LineWidth',3)

sigma = noise_level;
save('profile_data/BA2_noisy_data_2200.mat','t','output','ka1','ka2','kd1','kd2','R0s','AL20','Rmaxs','concs','t_asc','t_dis','sigma')

%%
kd2 = kd2_2;


for i=1:length(concs)
    t_star = t_stars(i);
    
    
    R0 = R0s(i);
    Am = concs(i);
    
    t_0 = t_asc(1) - t_star;
    t_pred_asc = t_0:1:(t_asc(1)-1);
    t_asc_temp = [t_pred_asc, t_asc];
    t_asc_temp = t_asc_temp - t_0;

    y0 = [Rmaxs(i), 0, 0];
    params = [ka1, ka2, kd1, kd2, Am];
    [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,params), t_asc_temp, y0, opts);
    y_asc(1:length(t_pred_asc),:) = [];
    

    y0 = y_asc(end,:);
    params = [ka1, ka2, kd1, kd2, 0];
    [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,params), t_dis, y0, opts);
    
    y_all(:,:,i) = [y_asc; y_dis(2:end,:)];

end

output = zeros(length(t),5);

for i=1:length(concs)
    output(:,i) = y_all(:,2,i) + y_all(:,3,i);
end

max_output = max(max(output));
for i=1:length(concs)
    output(:,i) = output(:,i) + normrnd(0,noise_level,[length(t),1]);
    scatter(t, output(:,i),10,'filled', 'MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410])
end
xline(t(151),'k--','LineWidth',2)
legend(strcat('k_{d2}=',num2str(kd2_1,'%.2e')),'','','','',strcat('k_{d2}=',num2str(kd2_2,'%.2e')),'','','','','','Location','Southeast')
[~,~] = title('Bivalent Analyte Model-2: Noisy Simulations for Extended Dissociation Length','FontSize',14);%;,...
%         {strcat(strcat('k_{a1}=',num2str(ka1,'%.2e')), strcat('-----k_{d1}=',num2str(kd1,'%.2e'))),...
%          strcat(strcat('k_{a2}=',num2str(ka2,'%.2e')), strcat('-----k_{d2}=',num2str(kd2,'%.2e'))),...
%          strcat(strcat('total\_t_{asc}=',num2str(t_asc(end) - t_asc(1))),      strcat('-----total\_t_{dis}=',num2str(t_dis(end) - t_dis(1))))
%          },
%      'FontSize',17);
set(gca,'FontSize',12)
box on
hold off

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
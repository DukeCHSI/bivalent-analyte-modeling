clear all; close all; clc;

table = readtable('../Bivalent Data/biphasic example newformat.xlsx');

t(:,1) = table.Ligand1_Concentration1;
y(:,1) = table.Var2;
t(:,2) = table.Ligand1_Concentration2;
y(:,2) = table.Var4;
t(:,3) = table.Ligand1_Concentration3;
y(:,3) = table.Var6;
t(:,4) = table.Ligand1_Concentration4;
y(:,4) = table.Var8;
t(:,5) = table.Ligand1_Concentration5;
y(:,5) = table.Var10;

%% Initialization
% Time interval [1,600]
t = 1:1:600;
% Initial conditions
y0 = [0.5; 0.1];
% Chosen parameters
params = [1e4, 1e-8, 1e-4];

% Pertubed percentage (1%)
h = 0.01;

% Intialize the pertubed parameters
params1 = params;
params2 = params;

% Get the number of parameters and the number of data
pars_change = [1, 2, 3];
num_par = length(pars_change);
num_ins = length(t);

% Initialize the two sensitivity matrix
% sum_S: sensitivity matrix for X1 + X2
% sep_S: sensitivity matrix for [X1, X2]
sum_S = zeros(num_ins, num_par);
sep_S = zeros(num_ins*length(y0), num_par);

%% Compute the sensitivity matrices

for ipar=1:num_par
    i = pars_change(ipar);
    % Pertube each parameter by 1% above and 1% below
    params1(i) = (1 + h)*params(i);
    params2(i) = (1 - h)*params(i);
    
    % Solve the ode with the new parameters
    [~,out1] = ode45(@(t,y) bivalent_rhs(t,y,params1), t, y0);
    [~,out2] = ode45(@(t,y) bivalent_rhs(t,y,params2), t, y0);
    
    % Calculate the change of sum using central finite difference
    sum_out1 = sum(out1,2);
    sum_out2 = sum(out2,2);
    sum_S(:,ipar) = (sum_out1-sum_out2)/(params1(i) - params2(i));
    
    % Calculate the change of [X1, X2] using central finite difference
    sep_out1 = [out1(:,1); out1(:,2)];
    sep_out2 = [out2(:,1); out2(:,2)];
    sep_S(:,ipar) = (sep_out1-sep_out2)/(params1(i) - params2(i));
    
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
    sum_cv(i) = sum_SE(i)/abs(params(i));
end

% R = mvnrnd(params,sum_C,1000);
% figure;
% scatter(R(:,1),R(:,2))
% xlabel('ka1')
% ylabel('kd1')
% figure;
% scatter(R(:,1),R(:,3))
% xlabel('ka1')
% ylabel('ka2')
% figure;
% scatter(R(:,1),R(:,4))
% xlabel('ka1')
% ylabel('kd2')
% figure;
% scatter(R(:,2),R(:,3))
% xlabel('kd1')
% ylabel('ka2')
% figure;
% scatter(R(:,2),R(:,4))
% xlabel('kd1')
% ylabel('kd2')
% figure;
% scatter(R(:,3),R(:,4))
% xlabel('ka2')
% ylabel('kd2')

%% Compute the number of identifiable parameters for [X1, X2]
% Remove the outputs at t=1
sep_S([1 1+num_ins],:) = [];

% Repeat the same steps for [X1, X2]
sep_FIM = sep_S'*sep_S;
sep_C = inv(sep_FIM);
invAA = sep_C*sep_FIM
sep_rank = rank(sep_FIM);

sep_SE = sqrt(diag(sep_C));

sep_cv = zeros(1,num_par);
for i=1:num_par
    sep_cv(i) = sep_SE(i)/abs(params(i));
end

%% bivalent right hand side model
function dy = bivalent_rhs(t,y, params)
X1 = y(1);
X2 = y(2);

B0 = 1.36;
C = 1e-7;

ka1 = params(1);
kd1 = 1e-4;
ka2 = params(2);
kd2 = params(3);

% ODE equations
dX2 = ka2*X1*(B0 - X1 - 2*X2) - 2*kd2*X2;
dX1 = 2*ka1*C*(B0 - X1 - 2*X2) - kd1*X1 - dX2;

dy = [dX1; dX2];
end
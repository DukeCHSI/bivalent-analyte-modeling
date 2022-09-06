clear all; close all; clc;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

%% Load data
T = readtable('test_data.csv');
T_asc_idx = T.AssocIndicator == 1;
T_dis_idx = T.DissocIndicator == 1;
data = T(T_asc_idx | T_dis_idx,:);
concs = unique(data.Concentration);



%%
load("best_results_data.mat")

params = best_params;

output = run_model_data(params,data, concs);




figure(2)
scatter(data.Time, data.RU)
hold on
scatter(output.Time,output.RU)
xlabel("Time")
ylabel("RU")
title("Attempt with 1 real data set")
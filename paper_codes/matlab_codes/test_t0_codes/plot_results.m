clear all; close all; clc;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
load('noisy_data.mat');

t_stars = [30, 50, 60, 70, 75];


params = [ka1, ka2, kd1, kd2, Rmaxs, t_stars];



true_y = run_model2(params,output,t,t_asc,t_dis, concs);

true_error = (true_y - output).^2;
true_error = sum(sum(true_error));

figure(1)
scatter(t,output)
hold on
plot(t,true_y)

%%
load("best_results_random_t0.mat")

params = best_params;

estimated_y = run_model2(params,output,t,t_asc,t_dis, concs);

estimated_error = (estimated_y - output).^2;
estimated_error = sum(sum(estimated_error));


figure(2)
scatter(t,output)
hold on
plot(t,estimated_y,'LineWidth',4)
xlabel('Time')
ylabel('RU')
title("Attempt with simulated data")

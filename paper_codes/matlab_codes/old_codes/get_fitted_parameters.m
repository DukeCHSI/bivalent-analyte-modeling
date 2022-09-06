clear all; close all; clc;

numT = 2200;

% save_dir = strcat('generated_data_N5_',num2str(numT),'T_monofirst/');
save_dir = 'generated_data_noise_5/';
load(strcat(save_dir,'data_5.mat'));

numIC = 81;


fitted_params_list = [];
IC_params_list = [];
fval_list = [];
for iIC = 1:numIC
    load(strcat(save_dir,'result_5/IC_',num2str(iIC),'/result.mat'));
    fitted_params_list = [fitted_params_list; params];
    IC_params_list = [IC_params_list; IC_params];
    fval_list = [fval_list; fval];
end
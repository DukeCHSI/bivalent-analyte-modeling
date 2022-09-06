clear all; close all; clc;

data = [];

load('generated_data_N5_1370T_oversampled/data_5.mat');


for i=1:length(concs)
    data = [data, [t', output(:,i)]];
end
%%
data = [];

load('generated_data_N5_1370T/data_5.mat');

for i=1:length(concs)
    data = [data, [t', output(:,i)]];
end

%%
data = [];
load('generated_data_N5_2200T_oversampled/data_5.mat');

for i=1:length(concs)
    data = [data, [t', output(:,i)]];
end

%%
data = [];
load('generated_data_N5_2200T/data_5.mat');

for i=1:length(concs)
    data = [data, [t', output(:,i)]];
end
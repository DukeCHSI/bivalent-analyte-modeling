clear all; close all; clc;

T = readtable('test_data.csv');

T_asc_idx = T.AssocIndicator == 1;
T_dis_idx = T.DissocIndicator == 1;

T_keep_idx = T_asc_idx + T_dis_idx;

T_asc = T(T_asc_idx,:);
T_dis = T(T_dis_idx,:);

concs = unique(T.Concentration);

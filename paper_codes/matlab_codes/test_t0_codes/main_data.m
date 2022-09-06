function main_data
opts = optimset('Display','iter','MaxIter',1e4,'TolFun',1e-8,'TolX',1e-8);

T = readtable('test_data.csv');

T_asc_idx = T.AssocIndicator == 1;
T_dis_idx = T.DissocIndicator == 1;

T_keep_idx = T_asc_idx + T_dis_idx;

data = T(T_asc_idx | T_dis_idx,:);

concs = unique(data.Concentration);

% Chosen parameters
lower_power = [1, -7, -7, -7];
upper_power = [6, -1, -1, -1];
lower_time = [0, 0, 0, 0, 0];
upper_time = [120, 120, 120, 120, 120];

num_run = 20;

best_fval = 1e64;

estimated_params = [];

Rmaxs = ones(1,length(concs))*max(data.RU);

for irun = 1:num_run
    k0 = 10.^(rand(1,4).*(upper_power - lower_power) + lower_power);
    tstar0 = rand(1,5).*(upper_time - lower_time) + lower_time;

    IC_params = [k0, Rmaxs, tstar0];
%     IC_params = [ka1, ka2, kd1, kd2, Rmaxs, 30, 50, 60, 70, 75];
    
    LB = zeros(size(IC_params),'like', IC_params);

    [opt_params,fval] = fmincon(@(params) objective_data(params,data,concs), IC_params,[],[],[],[],LB,[],[], opts);
    
    estimated_params = [estimated_params; IC_params, opt_params, fval];

    if fval < best_fval
        best_fval = fval;
        best_params = opt_params;
        best_inits = IC_params;
    end
end

save('best_results_data.mat','best_params','best_fval','best_inits','estimated_params')

end

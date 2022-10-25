function run_PI_new8
opts = optimset('Display','iter','MaxIter',1e4,'TolFun',1e-8,'TolX',1e-8);
rng(42);

concs = 2*[3.13e-08, 6.25e-08, 1.25e-07, 2.50e-07, 5.00e-07];
baseline = 120;
association = 300;
dissociation = 600;

num_run = 100;
num_data = 5;
data_name = 'data_new/dataset';
param_name = 'fitting_results_new_short/parameter';
save_name = 'profile_results_new_short/SSE';
N_sample = 20;
% Chosen lower/upper sampled values
lower_sample = [1e3, 1e-5, 1e-4, 1e-12];
upper_sample = [1e4, 1e-3, 1e-2, 5e-4];
% Generate the sampled values for the parameters
nump = length(lower_sample);
sampled_params = zeros(nump, N_sample);
for i=1:nump
    sampled_params(i,:) = 10.^(linspace(log10(lower_sample(i)),log10(upper_sample(i)),N_sample));
end

for idata = 8%:num_data

    dataset = load(strcat(data_name,num2str(idata),".mat"));
    dataset = struct2cell(dataset);
    dataset = dataset{1};
    Time = dataset{1};
    RU = dataset{2};
    Time = table2array(Time);
    RU = table2array(RU);

    temp_RU = zeros(size(RU));
    temp_Time = zeros(size(Time));
    for i=1:length(concs)
        temp = [Time(:,i), RU(:,i)];
        temp(any(isnan(temp), 2), :) = [];
        rows_to_remove = temp(:,1) < baseline;
        temp(rows_to_remove,:) = [];

        rows_to_remove = temp(:,1) > (baseline + association + dissociation);
        temp(rows_to_remove,:) = [];

        temp_Time(1:length(temp(:,1)),i) = temp(:,1);
        temp_RU(1:length(temp(:,2)),i) = temp(:,2);
    end

    RU = temp_RU;
    Time = temp_Time;

    clear temp_RU temp_Time

    estimated_params = load(strcat(param_name,num2str(idata),".mat"));
    full_SSE = estimated_params.min_SSE;
    full_params = estimated_params.best_params;

    clear estimated_params

    RUmaxs = max(max(RU))*ones(1,length(concs));
    LB = zeros(1,14);
    LB1 = zeros(1,14);

    sampled_params = [sampled_params'; full_params(1:4)]';
    sampled_params = sort(sampled_params,2);

    % Chosen parameters
    [~,N_sample] = size(sampled_params);
    SSE = zeros(nump,N_sample);
    for ip=1:nump
        for is=1:N_sample
            kinetics_pars = sampled_params(:,is);
            fixed_par = kinetics_pars(ip);
            % Remove fixed par
            kinetics_pars(ip) = [];
            LB(ip) = [];
            %
            min_SSE = 1e64;
            best_params = zeros(1,14);

            for ir=1:num_run
                formatSpec = 'IP is %f, IS is %f, IR is %f\n';
                fprintf(formatSpec,ip,is,ir)

                sampled_Rmaxs = RUmaxs.*rand(1,length(concs));
                sampled_tstars = 120.*rand(1,length(concs));

                IC_params = [kinetics_pars', sampled_Rmaxs, sampled_tstars];

                try
                    [min_params,min_fval] = fmincon(@(params) objective_PI(Time, RU, params, fixed_par, ip, ...
                        concs,baseline,association,dissociation), ...
                        IC_params,[],[],[],[],LB,[],[]);
                catch
                    min_fval = 4e64;
                    min_params = IC_params;
                end

                if min_fval < min_SSE
                    min_SSE = min_fval;
                    best_params = min_params;
                end
            end
            
            SSE(ip,is) = min_SSE;
            LB = LB1;
        end
    end
    save(strcat(save_name,num2str(idata),'.mat'),'sampled_params','SSE','full_SSE','full_params');
end

end

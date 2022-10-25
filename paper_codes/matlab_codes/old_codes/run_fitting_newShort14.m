function run_fitting_newShort14
opts = optimset('Display','iter','MaxIter',1e4,'TolFun',1e-8,'TolX',1e-8);
rng(42);

concs = 2*[3.13e-08, 6.25e-08, 1.25e-07, 2.50e-07, 5.00e-07];
baseline = 120;
association = 300;
dissociation = 600;

num_data = 14;
data_name = 'data_new/dataset';
save_name = 'fitting_results_short2/parameter';
for idata = 14%:num_data

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

    RUmaxs = max(max(RU))*ones(1,length(concs));
    LB = zeros(1,14);

    % Chosen parameters
    ka1_power = [2, 3, 4];
    ka2_power = [-6, -5, -4];
    kd1_power = [-4, -3, -2];
    kd2_power = [-6, -5, -4];
    sampled_params = [];
    for i1=1:3
        for i2=1:3
            for i3=1:3
                for i4=1:3
                    sampled_params = [sampled_params; ka1_power(i1), ka2_power(i2), kd1_power(i3), kd2_power(i4)];
                end
            end
        end
    end
    sampled_params = 10.^sampled_params;

    [num_sampled,~] = size(sampled_params);

    min_SSE = 1e64;
    best_params = zeros(1,14);
    for i=1:num_sampled
        i
        sampled_Rmaxs = RUmaxs.*rand(1,length(concs));
        sampled_tstars = 120.*rand(1,length(concs));

        IC_params = [sampled_params(i,:), sampled_Rmaxs, sampled_tstars];

        try
            [min_params,min_fval] = fmincon(@(params) objective(Time,RU,params,concs,baseline,association,dissociation), IC_params,[],[],[],[],LB,[],[]);
        catch
            min_fval = 4e64;
            min_params = IC_params;
        end
        
        if min_fval < min_SSE
            min_SSE = min_fval;
            best_params = min_params;
        end
    end

    save(strcat(save_name,num2str(idata),'.mat'),'best_params','min_SSE');
end

end
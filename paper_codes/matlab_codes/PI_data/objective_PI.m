function error = objective_PI(t,Ydata,params,fixed_par,ip,concs,baseline,association,dissociation)
params = insert(params, fixed_par, ip);
ka1 = params(1);
ka2 = params(2);
kd1 = params(3);
kd2 = params(4);
num_pars = 4;


error = 0;
num_data = 0;
for i=1:length(concs)
    data = [t(:,i), Ydata(:,i)];
    rows_to_remove = data(:,1) == 0;
    data(rows_to_remove,:) = [];
    rows_asc = data(:,1) < baseline + association;
    data_asc = data(rows_asc,:);
    data(rows_asc,:) = [];

    t_asc = data_asc(:,1);
    t_dis = data(:,1);
    Am = concs(i);

    Rmax = params(i+num_pars);
    t_star = params(i+5+num_pars);

    ode_params = [ka1, ka2, kd1, kd2, Am];
    y0 = [Rmax, 0, 0];

    if t_star >= 1
        t_0 = t_asc(1) - t_star;
        t_pred_asc = t_0:1:(t_asc(1)-1);
        t_asc_temp = [t_pred_asc'; t_asc];
        t_asc_temp = t_asc_temp - t_0;

        [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_asc_temp, y0);
        y_asc(1:length(t_pred_asc),:) = [];
    else
        [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_asc, y0);
    end

    y0 = y_asc(end,:);
    ode_params = [0, 0, kd1, kd2, 0];
    [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_dis, y0);

    y_pred = [y_asc(:,2) + y_asc(:,3); y_dis(:,2) + y_dis(:,3)];
    y_data = [data_asc(:,2); data(:,2)];
    error = error + sum((y_pred - y_data).^2);

    if length(y_data) ~= length(y_pred)
        disp("Lengths are not equal")
    end

    num_data = num_data + length(y_data);
end
error = error;%/num_data;
end

function err = objective(params,data,t,t_asc,t_dis, concs)
ode_params1 = params(1:4);
ode_params = ode_params1;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

y = zeros(length(t),length(concs));



% Solve the ode
for i=1:length(concs)
    Am = concs(i);
    Rmax = params(4+i);
    ode_params(5) = Am;
    
    t_star = params(4+length(concs)+i);

    if t_star >= 1
        t_0 = t_asc(1) - t_star;
        t_pred_asc = t_0:1:(t_asc(1)-1);
        t_asc_temp = [t_pred_asc, t_asc];
        t_asc_temp = t_asc_temp - t_0;
        y0 = [Rmax, 0, 0];
        [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_asc_temp, y0, opts);
        y_asc(1:length(t_pred_asc),:) = [];
    else
        y0 = [Rmax, 0, 0];
        [~, y_asc] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_asc, y0, opts);
    end

    y0 = y_asc(end,:);
    ode_params(5) = 0;
    [~, y_dis] = ode15s(@(t,y) bivalent_rhs(t,y,ode_params), t_dis, y0, opts);
    
    y(:,i) = [y_asc(:,2) + y_asc(:,3); y_dis(2:end,2) + y_dis(2:end,3)];
    
    ode_params = ode_params1;
end

err = sum(sum((y-data).^2));
end

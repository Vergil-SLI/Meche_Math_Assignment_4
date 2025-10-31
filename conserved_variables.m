function conserved_variables()
    % Preset values / values to keep consistent across comparison
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;    
    V0 = [x0;y0;dxdt0;dydt0];
    h_ref = 0.05;
    rate_func_in = @(t, V) gravity_rate_func(t,V, orbit_params);

    
    % heun's method struct values__________________________________________
    BT_struct = struct();
    BT_struct.A = [0, 0; 1, 0]; % matrix of a_{ij} values
    BT_struct.B = [0.5, 0.5];% vector of b_i values
    BT_struct.C = [0, 1]; % vector of c_i values
    [t_list,V_list,~, ~] = explicit_RK_fixed_step_integration(rate_func_in,[0,30],V0,h_ref,BT_struct);

    x_vals = V_list(1, :);
    y_vals = V_list(2, :);
    dx_vals = V_list(3, :);
    dy_vals = V_list(4, :);

    % calculate the mechanical energy, E, and angular momentum, H
    mag_r = sqrt(x_vals.^2 + y_vals.^2);
    E_heun = (1/2) * orbit_params.m_planet * (dx_vals.^2 + dy_vals.^2) - ...
        (orbit_params.m_sun * orbit_params.m_planet * orbit_params.G) ./ mag_r;
    H_heun = orbit_params.m_planet * (x_vals .* dy_vals - y_vals .* dx_vals);
    
    % explicit midpoint method struct values________________________________
    BT_struct = struct();
    BT_struct.A = [0, 0; 0.5, 0]; % matrix of a_{ij} values
    BT_struct.B = [0, 1];% vector of b_i values
    BT_struct.C = [0, 0.5]; % vector of c_i values
    [t_list,V_list,~, ~] = explicit_RK_fixed_step_integration(rate_func_in,[0,30],V0,h_ref,BT_struct);

    x_vals = V_list(1, :);
    y_vals = V_list(2, :);
    dx_vals = V_list(3, :);
    dy_vals = V_list(4, :);

    % calculate the mechanical energy, E, and angular momentum, H
    mag_r = sqrt(x_vals.^2 + y_vals.^2);
    E_midpoint = (1/2) * orbit_params.m_planet * (dx_vals.^2 + dy_vals.^2) - ...
        (orbit_params.m_sun * orbit_params.m_planet * orbit_params.G) ./ mag_r;
    H_midpoint = orbit_params.m_planet * (x_vals .* dy_vals - y_vals .* dx_vals); 
    
    % forward euler method struct values________________________________
    BT_struct = struct();
    BT_struct.A = [0]; % matrix of a_{ij} values
    BT_struct.B = [1];% vector of b_i values
    BT_struct.C = [0]; % vector of c_i values
    [t_list,V_list,~, ~] = explicit_RK_fixed_step_integration(rate_func_in,[0,30],V0,h_ref,BT_struct);

    x_vals = V_list(1, :);
    y_vals = V_list(2, :);
    dx_vals = V_list(3, :);
    dy_vals = V_list(4, :);

    % calculate the mechanical energy, E, and angular momentum, H
    mag_r = sqrt(x_vals.^2 + y_vals.^2);
    E_euler = (1/2) * orbit_params.m_planet * (dx_vals.^2 + dy_vals.^2) - ...
        (orbit_params.m_sun * orbit_params.m_planet * orbit_params.G) ./ mag_r;
    H_euler = orbit_params.m_planet * (x_vals .* dy_vals - y_vals .* dx_vals); 

    
    % plotting mechanical energy & angular momentum
    figure(1)
    plot(t_list, E_heun)
    title("mechanical energy vs. time (href = 0.05)")
    hold on
    plot(t_list, E_midpoint)
    plot(t_list, E_euler)
    legend("heun's method", "explicit midpoint", "forward euler")
    ylim([-4, -1])

    figure(2)
    plot(t_list, H_heun)
    title("angular momentum vs. time (href = 0.05)")
    hold on
    plot(t_list, H_midpoint)
    plot(t_list, H_euler)
    legend("heun's method", "explicit midpoint", "forward euler")
    ylim([11.5, 15.5])
end
%this function computes the orbit of a planet about a sun
%the sun is assumed to located at the origin (0,0)
%and motion is restricted to the x-y plane
%INPUTS:
%t_list: a list of times to compute the position & velocity of the planet
%V0: initial condition. V0 is a column vector consisting
% of the initial position and velocity of the planet:
% V0 = [x(0); y(0); dx/dt(0); dy/dt(0)]
%orbit_params: a struct describing the system parameters
% orbit_params.m_sun: mass of the sun
% orbit_params.m_planet: mass of the planet
% orbit_params.G: gravitational constant
% Force = -m_planet*m_sun*G/rˆ2
%OUTPUTS:
%V_list: the state of the planet at each time in t_list
% if t_list is a list, then V_list is a Nx4 MATRIX
% where each ROW has the form [x_i,y_i,dx/dt_i,dy/dt_i]
% corresponding to the values at the ith time
% if t_list is a SCALAR (i.e. t_list = t),
% then V_list is a COLUMN VECTOR of the form:
% [x(t); y(t); dx/dt(t); dy/dt(t)]
%NOTES:
%This function needs all the other functions in this file to run
%I HIGHLY RECOMMEND JUST SAVING THIS FUNCTION IN ITS OWN FILE
%DON'T COPY AND PASTE INTO YOUR CODE! IT'S NOT WORTH IT!
%
%USAGE EXAMPLE:
%At the bottom of this file is a function usage_example()
%which shows how to use compute_planetary_motion(...)
%You can start from there and then tweak it.

function usage_example()
    % Actual solution calc
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,30,100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);

    % plot actual solution
    axis equal; axis square;
    axis([-20,20,-20,20])
    hold off
    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    hold on
    plot(V_list(:,1),V_list(:,2),'k');

    % RK approximation calculations________________________________________
    % % forward euler
    % h_ref = 0.01;
    % BT_struct = struct();
    % BT_struct.A = [0]; % matrix of a_{ij} values
    % BT_struct.B = [1];% vector of b_i values
    % BT_struct.C = [0]; % vector of c_i values
    % 
    % rate_func_in = @(t, V) gravity_rate_func(t,V, orbit_params);
    % 
    % [t_list,V_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in,[0,30],V0,h_ref,BT_struct);
    % x_vals = V_list(1, :);
    % y_vals = V_list(2, :);
    % 
    % 
    % plot(x_vals, y_vals, "g")
    % title("Comparing Forward Euler Approx. to True Solution (h = 0.01)")
    % legend("", "true solution", "forward euler")

    % % explicit midpoint
    h_ref = 0.01;
    BT_struct = struct();
    BT_struct.A = [0, 0; 0.5, 0]; % matrix of a_{ij} values
    BT_struct.B = [0, 1];% vector of b_i values
    BT_struct.C = [0, 0.5]; % vector of c_i values

    % DormandPrince = struct();
    % DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    % DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
    % 5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    % DormandPrince.A = [0,0,0,0,0,0,0;
    %     1/5, 0, 0, 0,0,0,0;...
    %     3/40, 9/40, 0, 0, 0, 0,0;...
    %     44/45, -56/15, 32/9, 0, 0, 0,0;...
    %     19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
    %     9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
    %     35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];
    % 
    % BT_struct = DormandPrince;

    rate_func_in = @(t, V) gravity_rate_func(t,V, orbit_params);

    % [t_list,V_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in,[0,30],V0,h_ref,BT_struct);
    % x_vals = V_list(1, :);
    % y_vals = V_list(2, :);
    % 
    % 
    % plot(x_vals, y_vals, "r")
    % title("Comparing Explicit Midpoint Approx. to True Solution (h = 0.01)")
    % legend("", "true solution", "explicit midpoint")
    
    
    % % heun's method
    % h_ref = 0.01;
    % BT_struct = struct();
    % BT_struct.A = [0, 0; 1, 0]; % matrix of a_{ij} values
    % BT_struct.B = [0.5, 0.5];% vector of b_i values
    % BT_struct.C = [0, 1]; % vector of c_i values
    % 
    % rate_func_in = @(t, V) gravity_rate_func(t,V, orbit_params);

    [t_list,V_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in,[0,30],V0,h_ref,BT_struct);
    x_vals = V_list(1, :);
    y_vals = V_list(2, :);
    dx_vals = V_list(3, :);
    dy_vals = V_list(4, :);


    plot(x_vals, y_vals, "b")
    title("Comparing Heun's Approx. to True Solution (h = 0.01)")
    legend("", "true solution", "heun's method")
    

    mag_r = sqrt(x_vals.^2 + y_vals.^2);
    E = (1/2) * orbit_params.m_planet * (dx_vals.^2 + dy_vals.^2) - ...
        (orbit_params.m_sun * orbit_params.m_planet * orbit_params.G) ./ mag_r;
    H = orbit_params.m_planet * (x_vals .* dy_vals - y_vals .* dx_vals);

    hold off
    plot(t_list, E)
    hold on
    plot(t_list, H)
    legend("mechanical energy", "angular momentum")
    hold off

    n_samples = 30;
    h_ref_list = logspace(-3,1, 30);
    
    num_evals_list = zeros(1, n_samples);
    h_avg_list = zeros(1, n_samples);
    tr_error_list = zeros(1, n_samples);
    V_list = compute_planetary_motion(t_range,V0,orbit_params)';

    % Local truncation error of Euler and Explicit Midpoint Method
    t_ref = 0.1;
    hspan = [-5, -1, 100];
    [h_list, analytical_difference,expmid_error_list] = local_truncation_error_assignment_4(t_ref, hspan, rate_func_in);
    
    % Find line of best fit coefficients and estimate p-value
    filter_params = struct();
    filter_params.min_yval = 1e-14;
    filter_params.max_yval = 1e-5;
   
    [ana_p, ~] = loglog_fit(h_list, analytical_difference);
    fprintf("Local analytical p-value: ");
    disp(ana_p);
    [expmid_p,expmid_k] = loglog_fit(h_list,expmid_error_list,filter_params);
    expmid_y_data = expmid_k.*((h_list.^expmid_p));
    fprintf("Local explicit Midpoint p-value: ");
    disp(expmid_p);

    % Plot log scale graph of errors vs h_list (all the different step sizes used)
    figure(1);
    loglog(h_list,expmid_error_list,'go','MarkerFaceColor','g'); hold on
    loglog(h_list,expmid_y_data,'k','LineWidth',2); hold on
    loglog(h_list, analytical_difference, 'r-'); hold on

    % Set axes and legend
    xlim([1e-5 1e0])
    title("Local truncation error of explicit midpoint method vs. h list");
    xlabel("Time step sizes");
    ylabel("Error");
    legend("Explicit Midpoint Method", "Fit line",  "Analytical Difference");

    for n = 1:length(h_ref_list)
        h_ref = h_ref_list(n);

        [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in,[0,30],V0,h_ref,BT_struct);

        tr_error = norm(X_list(:,end) - V_list(:,end));
        tr_error_list(n) = tr_error;
        h_avg_list(n) = h_avg;
        num_evals_list(n) = num_evals;
    end

    filter_params = struct();
    filter_params.min_yval = 1e-10;
    filter_params.max_yval = 1e-2;
    
    [p1,k1] = loglog_fit(h_avg_list,tr_error_list, filter_params);
    [p2,k2] = loglog_fit(num_evals_list, tr_error_list, filter_params);

    
    p1
    p2
    figure(2)
    loglog(h_avg_list, tr_error_list, 'ro', 'markerfacecolor', 'r');
    hold on;
    loglog(h_avg_list, k1*h_avg_list.^p1, 'r-', 'markerfacecolor', 'r');

    figure(3)
    loglog(num_evals_list, tr_error_list, 'bo', 'markerfacecolor', 'b');
    hold on
    loglog(num_evals_list, k2*num_evals_list.^p2, 'b-', 'markerfacecolor', 'r');

end
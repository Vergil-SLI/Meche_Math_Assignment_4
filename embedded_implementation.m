function embedded_implementation()
    %Local truncation error experiments for embedded method
    % clear all
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
    V0 = [x0;y0;dxdt0;dydt0];
    n_samples = 60;
    h_ref_list = logspace(-3.3, 1, n_samples);
    tspan = [0];
    abs_diff_list = zeros(1, n_samples);
    xb_diff_list = zeros(1, n_samples);
    tr_error_list1 = zeros(1, n_samples);
    tr_error_list2 = zeros(1, n_samples);

    rate_func_in = @(t, V) gravity_rate_func(t,V, orbit_params);
    
    % set up embedded method struct
    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
    5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
        1/5, 0, 0, 0,0,0,0;...
        3/40, 9/40, 0, 0, 0, 0,0;...
        44/45, -56/15, 32/9, 0, 0, 0,0;...
        19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
        9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
        35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

    BT_struct = DormandPrince;

    for n = 1:length(h_ref_list)
        h_ref = h_ref_list(n);
        V_list = compute_planetary_motion(tspan(1)+h_ref, V0, orbit_params);
    
        [XB1, XB2, ~] = RK_step_embedded(rate_func_in, tspan(1)+h_ref, V0, h_ref, BT_struct);
        abs_diff_list(n) = norm(V_list-V0);
        xb_diff_list(n) = norm(XB1-XB2);
        tr_error_list1(n) = norm(XB1 - V_list);
        tr_error_list2(n) = norm(XB2 - V_list);
    
    end
    
    filter_params = struct();
    filter_params.min_yval = 1e-13;
    filter_params.max_yval = 1e-6;
    
    [p1, k1] = loglog_fit(h_ref_list, tr_error_list1, filter_params);
    [p2, k2] = loglog_fit(h_ref_list, tr_error_list2, filter_params);
    
    p1
    p2
    

    % hold off
    fig1 = figure(3)
    set(0,'CurrentFigure',fig1)
    loglog(h_ref_list, tr_error_list1, 'bo', 'markerfacecolor', 'b', 'markersize', 3);
    hold on
    loglog(h_ref_list, tr_error_list2, 'ro', 'markerfacecolor', 'r', 'markersize', 3);
    loglog(h_ref_list, xb_diff_list, 'go', 'markerfacecolor', 'g', 'markersize', 3);
    loglog(h_ref_list, abs_diff_list, 'mo', 'markerfacecolor', 'm', 'markersize', 3);
    xlabel("Timestep h")
    lgd = legend("XB1 Local Truncation Error", "XB2 Local Truncation Error", "|XB1 - XB2|", "|f(tref + h) - f(t)}");
    lgd.Location = 'southeast';

    [p1, k1] = loglog_fit(h_ref_list, tr_error_list1, filter_params);
    [p2, k2] = loglog_fit(h_ref_list, tr_error_list2, filter_params);
    [p3, k3] = loglog_fit(h_ref_list, xb_diff_list, filter_params);
    
    p1
    p2
    p3

    % hold off
    fig2 = figure(2);
    set(0,'CurrentFigure',fig2)
    loglog(xb_diff_list, tr_error_list1, 'bo', 'markerfacecolor', 'b', 'markersize', 3);
    hold on
    loglog(xb_diff_list, tr_error_list2, 'ro', 'markerfacecolor', 'r', 'markersize', 3);
    lgd = legend("XB1 Local Truncation Error", "XB2 Local Truncation Error");
    xlabel("|XB1-XB2|")
    lgd.Location = 'southeast';

    [p4, k4] = loglog_fit(xb_diff_list, tr_error_list1, filter_params);
    [p5, k5] = loglog_fit(xb_diff_list, tr_error_list2, filter_params);

    p4
    p5
end
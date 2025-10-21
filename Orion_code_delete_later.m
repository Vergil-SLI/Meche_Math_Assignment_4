%Local truncation error experiments for embedded method
n_samples = 60;
h_ref_list = logspace(-6, 1, n_samples);

abs_diff_list = zeros(1, n_samples);
tr_error_list1 = zeros(1, n_samples);
tr_error_list2 = zeros(1, n_samples);

for n = 1:length(h_ref_list)
    h_ref = h_ref_list(n);
    V_list = compute_planetary_motion(tspan(1)+h_ref, V0, orbit_params);

    [XB1, XB2, ~] = explicit_RK_step_embedded(my_rate, tspan(1), V0, h_ref, BT_struct);
    abs_diff_list(n) = norm(V_list-V0);
    tr_error_list1(n) = norm(XB1 - V_list);
    tr_error_list2(n) = norm(XB2 - V_list);

end

filter_params = struct();
filter_params.min_yval = 1e-13;
filter_params.max_yval = 1e-6;

[p1, k1] = loglog_fit(h_ref_list, tr_error_list1, filterparams);
[p2, k2] = loglog_fit(h_ref_list, tr_error_list2, filterparams);

p1
p2

figure(2);
loglog(h_ref_list, abs_diff_list, 'ro', 'markerfacecolor', 'r', 'markersize', 3);
hold on
loglog(h_ref_list, tr_error_list1, 'bo', 'markerfacecolor', 'r', 'markersize', 3);
loglog(h_ref_list, tr_error_list2, 'go', 'markerfacecolor', 'r', 'markersize', 3);


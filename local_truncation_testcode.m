% local truncation error experiments for embedded methods
n_samples = 60;
h_ref_list = logspace(-6,1, n_samples)
tr_error_list1 = [];
tr_error_list2 = [];

for n = 1:length(h_ref_list)
    h_ref = h_ref_list(n)
    V_list = compute_planetary_motion(t_range,V0,orbit_params);

    [XB1, XB2, ~] = explicit_RK_step_embedded(rate_func, tspan(1), V0, h_ref, DormandPrince)

    abs_diff_list(n) = norm(V_list-V0);
    tr_error_list1(n) = norm(XB1-V_list);
    tr_error_list2(n) = norm(XB2-V_list);
end

filter_params = struct();
filter_params.min_yval = 1e-13;
filter_params.max_yval = 1e-6;
    
[p1,k1] = loglogfit(h_ref_list,tr_error_list1, filter_params);
[p2,k2] = loglogfit(h_ref_list, tr_error_list2, filter_params);

figure();
loglog(h_ref_list, abs_diff_list, 'ro', 'markerfacecolor', 'r');
hold on
loglog(h_ref_list, tr_error_list1, 'bo', 'markerfacecolor', 'b');
loglog(h_ref_list, tr_error_list2, 'go', 'markerfacecolor', 'g');
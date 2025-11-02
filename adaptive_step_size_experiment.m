% starting values
clear all
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;
x0 = 10;
y0 = 5;
dxdt0 = 0;
dydt0 = 1.5;
V0 = [x0;y0;dxdt0;dydt0];
rate_func_in = @(t, V) gravity_rate_func(t,V, orbit_params);
tspan = [0,30];
p = 5;

t_range = linspace(0,30,100);
V_list = compute_planetary_motion(t_range,V0,orbit_params)';

% adaptive step size
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
h_ref = 0.05;
error_desired_list = logspace(-15,0, 30);

for n = 1:length(error_desired_list)
    error_desired = error_desired_list(n);
    [t_list,X_list,h_avg, step_failure_rate, num_evals] = explicit_RK_variable_step_integration(rate_func_in,tspan,V0,h_ref,DormandPrince,p,error_desired);

    adaptive_tr_error = norm(X_list(:,end) - V_list(:,end));
    adaptive_tr_error_list(n) = adaptive_tr_error;
    adaptive_h_avg_list(n) = h_avg;
    adaptive_num_evals_list(n) = num_evals;
    adaptive_step_failure_rate(n) = step_failure_rate;
end


% fixed step size_________________________________________________________
BT_struct = struct();
BT_struct.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
BT_struct.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
BT_struct.A = [0,0,0,0,0,0,0;
    1/5, 0, 0, 0,0,0,0;...
    3/40, 9/40, 0, 0, 0, 0,0;...
    44/45, -56/15, 32/9, 0, 0, 0,0;...
    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
    9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

h_ref_list = logspace(-3,1, 30);
for n = 1:length(h_ref_list)
    h_ref = h_ref_list(n);

    [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in,tspan,V0,h_ref,BT_struct);

    fixed_tr_error = norm(X_list(:,end) - V_list(:,end));
    fixed_tr_error_list(n) = fixed_tr_error;
    fixed_h_avg_list(n) = h_avg;
    fixed_num_evals_list(n) = num_evals;
end

% loglog plots
figure(1)
loglog(adaptive_tr_error_list, adaptive_h_avg_list, "bo", MarkerFaceColor="b")
hold on
loglog(fixed_tr_error_list, fixed_h_avg_list, "ro", MarkerFaceColor="r")
legend("adaptive step size", "fixed step size")
title("global truncation error vs. avg. step size")
xlabel("avg step size")
ylabel("global truncation error")

figure(2)
loglog(adaptive_tr_error_list, adaptive_num_evals_list, "bo", MarkerFaceColor="b")
hold on
loglog(fixed_tr_error_list, fixed_num_evals_list, "ro", MarkerFaceColor="r")
legend("adaptive step size", "fixed step size")
title("global truncation error vs. # function evals")
xlabel("# of function evals")
ylabel("global truncation error")

% how often things fail
figure(3)
semilogx(adaptive_h_avg_list, adaptive_step_failure_rate, "bo", MarkerFaceColor="b")
title("failure rate as a function of avg. step size")
xlabel("avg step size")
ylabel("failure rate")

% plotting planetary path
error_desired = 0.0001;
[t_list,X_list,h_avg, step_failure_rate, num_evals] = explicit_RK_variable_step_integration(rate_func_in,tspan,V0,h_ref,DormandPrince,p,error_desired);
x_vals = X_list(1, :);
y_vals = X_list(2, :);
dx_vals = X_list(3, :);
dy_vals = X_list(4, :);
t_vals = t_list;

figure(4)
plot(x_vals, t_vals,'ro-','markerfacecolor','k','markeredgecolor','k','markersize',2)
hold on
plot(y_vals, t_vals, 'bo-','markerfacecolor','k','markeredgecolor','k','markersize',2)
title("position over time")
legend("x pos", "y pos")
xlabel("time")
ylabel("position")

figure(5)
plot(dx_vals, t_vals,'ro-','markerfacecolor','k','markeredgecolor','k','markersize',2)
hold on
plot(dy_vals, t_vals, 'bo-','markerfacecolor','k','markeredgecolor','k','markersize',2)
title("velocity over time")
legend("x vel", "y vel")
xlabel("time")
ylabel("velocity")

% whatever this scatter plot distance thing is
figure(6)
h_vals = diff(t_vals);
dist_between_planets = sqrt(x_vals.^2 + y_vals.^2);
plot(dist_between_planets(2:end), h_vals, "bo", MarkerFaceColor="b")
title("step size in relation to distance between the planet and the sun")
xlabel("distance")
ylabel("step size")
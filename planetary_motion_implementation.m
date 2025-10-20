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

    % RK approximation calculations
    h_ref = 0.01;
    BT_struct = struct();
    BT_struct.A = [0]; % matrix of a_{ij} values
    BT_struct.B = [1];% vector of b_i values
    BT_struct.C = [0]; % vector of c_i values
    
    rate_func_in = @(V) gravity_rate_func(0,V,orbit_params);
    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in,[0,30],V0,h_ref,BT_struct);

    
    % plot actual solution
    % axis equal; axis square;
    % axis([-20,20,-20,20])
    % hold on
    % plot(0,0,'ro','markerfacecolor','r','markersize',5);
    % plot(V_list(:,1),V_list(:,2),'k');
end
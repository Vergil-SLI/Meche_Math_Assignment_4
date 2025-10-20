%Rate function describing Newton's law of gravitation
%INPUTS:
%t: the current time
%V: the vector of the position and velocity of the planet
% V = [x_p; y_p; dxdt_p; dydt_p]
%orbit_params: a struct describing the system parameters
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant
%OUTPUTS:
%dVdt: a column vector describing the time derivative of V:
% dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]
function dVdt = gravity_rate_func(t,V,orbit_params)
    X = [V(1), V(2); V(3), V(4)];
    mag_r = sqrt(X(1,1)^2 + X(1,2)^2);
    ms = orbit_params.m_sun;
    mp = orbit_params.m_planet;
    G = orbit_params.G;

    k = -1* (mp*ms*G)/(mag_r^3);
    A = [0 1; k/mp, 0];

    temp = A * X;
    dVdt = [temp(1,1); temp(1,2); temp(2,1); temp(2,2)];
end


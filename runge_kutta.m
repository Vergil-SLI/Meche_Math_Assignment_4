%This function computes the value of X at the next time step
%for any arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB, num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct)
    k_vals = zeros(1, length(BT_struct.C));

    for i = 1:length(BT_struct.C)
        t_input = t + BT_struct.C(i) * h;
        X_input = 0;       
        for j = 1:i-1
            X_input = X_input + BT_struct.A(i, j) * k_vals(j);
        end
        X_input = XA + h*X_input;

        k_vals(i) = rate_func_in(t_input, X_input);
    end
end
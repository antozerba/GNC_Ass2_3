function [Ad, Bd, Cd, Ed] = discrete_system(K,T, h)
% DISCRETE_SYSTEM  Discrete-time model for the heading Kalman filter
%                  (Fossen 2021, Ch. 13.4), based on the first-order
%                  Nomoto model with a rudder bias modelled as a random walk.
%
%   [Ad, Bd, Cd, Ed] = discrete_system(K, T, h)
%
%   Continuous-time model:
%       psi_dot = r
%       r_dot   = -(1/T) r + (K/T)(delta - b) + w_r
%       b_dot   = w_b
%       y       = psi + v
%
%   INPUTS
%     K  - Nomoto gain [1/s]
%     T  - Nomoto time constant [s]
%         
%     h  - Sampling time [s].
%
%   OUTPUTS
%     Ad - Discrete state matrix (3x3).
%     Bd - Discrete input matrix (3x1).
%     Cd - Measurement matrix (1x3): only the yaw angle is measured.
%     Ed - Discrete process-noise matrix (3x2).
    A = [   
        0 1 0;
        0 -1/T -K/T;
        0 0 0;
        ];
    B = [0 K/T 0 ]';
    E = [
        0 0;
        1 0;
        0 1;
        ];
    C = [1 0 0];
    
    [Ad,Bd] = c2d(A,B,h);
    [Ad,Ed] = c2d(A,E,h);
    Cd = C;
end
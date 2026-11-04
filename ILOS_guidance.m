function [chi_d, y_int] = ILOS_guidance(e_y,pi_h, y_int, h)
% ILOS_GUIDANCE  Integral Line-of-Sight guidance law for straight-line path
%                following (Fossen 2021, Eqs. 12.108-12.109).
%
%   [chi_d, y_int] = ILOS_guidance(e_y, pi_h, y_int, h)

%   INPUTS
%     e_y    - Cross-track error [m]: lateral distance from the craft to
%              the active path segment, positive to starboard.
%     pi_h   - Path-tangential angle [rad]: direction of the active
%              segment from waypoint p_k to p_k+1, measured from North.
%     y_int  - Integral state [s] at the current time step. Reset to zero
%              by the caller at every waypoint switch.
%     h      - Sampling time [s], used for Euler integration of y_int.
%
%   OUTPUTS
%     chi_d  - Desired COURSE angle [rad].
%     y_int  - Integral state [s] updated to the next time step.


    k = 0.95;                   %desing param for integral gain
    L = 161;                    % ship length (m)
    delta_h = 5*L;              % lookhaed distance
    Kp = 1/delta_h;
    Ki = k*Kp;
    chi_d = pi_h - atan(Kp*e_y + Ki*y_int);

    y_int_dot = delta_h*e_y / (delta_h^2 + (e_y + k*y_int)^2);
    y_int = y_int + h * y_int_dot;



    
end
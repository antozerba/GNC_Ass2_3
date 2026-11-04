function [chi_d,last,e_y] = guidance(x,y,h)
% GUIDANCE  Waypoint-based ILOS guidance for straight-line path following.
%   y_int is a persistent integral state, reset at every waypoint switch.
%   NB: call "clear guidance" at the start of main.m.

    law = 'ilos';

    persistent y_int_p pi_h_prev
    if isempty(y_int_p)
        y_int_p   = 0;
        pi_h_prev = NaN;
    end
    
    [xk1,yk1,xk,yk,last] = WP_selector(x,y);
    [e_y,pi_h]           = crossTrackError(xk1,yk1,xk,yk,x,y);
    
    if pi_h ~= pi_h_prev          % in tnew segment reset integral state
        y_int_p = 0;
    end

    pi_h_prev = pi_h;
    
    if strcmp(law,'ilos')
        [chi_d, y_int_p] = ILOS_guidance(e_y, pi_h, y_int_p, h);
    else
        chi_d = LOS_guidance(e_y, pi_h);
    end
            
end
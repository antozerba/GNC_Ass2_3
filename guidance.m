function [chi_d, last] = guidance(x,y)
% GUIDANCE  Waypoint-based LOS guidance system for straight-line path
%           following (Fossen 2021, Ch. 12.4)
%   [chi_d, last] = guidance(x, y)
%
%   INPUTS
%     x, y   - Current North and East position 
%
%   OUTPUTS
%     chi_d  - Desired COURSE angle [rad].
%     last   - true when the craft is on the final path segment


    [xk1,yk1,xk,yk,last] = WP_selector(x,y);
    [e_y, pi_h] = crossTrackError(xk1,yk1,xk,yk,x,y); 
    chi_d = LOS_guidance(e_y, pi_h);          
    
end
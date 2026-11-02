function chi_d = LOS_guidance(e_y,pi_h)
% LOS_GUIDANCE  Lookahead-based Line-of-Sight guidance law for straight-line
%               path following (Fossen 2021, Eq. 12.78).
%   INPUTS
%     e_y   - Cross-track error [m]: lateral distance from the craft to the
%             active path segment, expressed in the path-fixed frame.
%      
%     pi_h  - Path-tangential angle [rad]: direction of the active segment
%             from waypoint p_k to p_k+1, measured from North
%
%   OUTPUT
%     chi_d - Desired COURSE angle [rad] 

    L = 161;                % ship length (m)
    delta_h = 5*L;          % lookhaed distance
    K = 1/delta_h;
    chi_d = pi_h - atan(K*e_y);
end
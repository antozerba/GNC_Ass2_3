function [e_y,pi_h] = crossTrackError(xk1,yk1,xk,yk,x,y)
% CROSSTRACKERROR  Computes the path-tangential angle of the active
%                  straight-line segment and the cross-track error of the
%                  craft relative to that segment (Fossen 2021, Ch. 12.4).
%   [e_y, pi_h] = crossTrackError(xk1, yk1, xk, yk, x, y)

%   INPUTS
%     xk1, yk1 - North and East coordinates [m] of the NEXT waypoint p_k+1
%                (end of the active segment).
%     xk,  yk  - North and East coordinates [m] of the PREVIOUS waypoint p_k
%                (start of the active segment).
%     x,   y   - Current North and East position of the craft [m] in NED.
%
%   OUTPUTS
%     e_y  - Cross-track error [m] distance from the craft to the (infinite) line through p_k and p_k+1.
%     pi_h - Path-tangential angle [rad] in (-pi, pi]: direction of the
%            segment from p_k to p_k+1


    pi_h = atan2(yk1-yk, xk1-xk);                    % path-tangential angle [rad]
    e_y  = -(x-xk)*sin(pi_h) + (y-yk)*cos(pi_h);     % cross-track error [m]
end
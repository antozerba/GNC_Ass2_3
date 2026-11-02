function [xk1,yk1,xk,yk,last] = WP_selector(x,y)
% WP_SELECTOR  Selects the active straight-line path segment p_k -> p_k+1
%              using a circle-of-acceptance criterion (Fossen 2021, Ch. 12.4).
%
%   [xk1, yk1, xk, yk, last] = WP_selector(x, y)
%
%
%   INPUTS
%     x, y  - Current North and East position of the craft [m] in NED.
%
%   OUTPUTS
%     xk,  yk   - North/East coordinates [m] of the previous waypoint p_k
%     xk1, yk1  - North/East coordinates [m] of the next waypoint p_k+1
%     last      - true when the active segment is the final one
%
%   PARAM
%     R  - Radius of the circle of acceptance [m]

persistent k WP R
if isempty(k)
    data = load('WP.mat');
    WP   = data.WP;        % 2 x n, row 1 = North, row 2 = East
    k    = 1;              % index of the active segment
    R    = 1000;           % radius of circle of acceptance [m]
end
n = size(WP,2);            % number of waypoints

% Active segment
xk  = WP(1,k);   yk  = WP(2,k);
xk1 = WP(1,k+1); yk1 = WP(2,k+1);

% Circle-of-acceptance test on the next waypoint
dx    = xk1 - x;
dy    = yk1 - y;
accpt = (dx^2 + dy^2 <= R^2);

% Switch to the next segment 
if accpt && k < n-1
    k   = k+1;
    xk  = WP(1,k);   yk  = WP(2,k);
    xk1 = WP(1,k+1); yk1 = WP(2,k+1);
end

last = (k == n-1);         %last seg
end
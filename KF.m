function [x_pst, P_pst, x_prd, P_prd, innov, S] = KF(x_prd, P_prd, Ad, Bd, Ed, Cd, Qd, Rd, psi_meas, delta)
%   [x_pst, P_pst, x_prd, P_prd, innov, S] = KF(x_prd, P_prd, Ad, Bd, Ed, Cd, Qd, Rd, psi_meas, delta)
%
%
%   INPUTS
%     x_prd    - Predicted state at step k (3x1): [psi; r; b] [rad, rad/s, rad].
%     P_prd    - Predicted error covariance at step k (3x3).
%     Ad,Bd,Ed - Discrete-time model matrices (from discrete_system).
%     Cd       - Measurement matrix (1x3), Cd = [1 0 0].
%     Qd       - Process noise covariance (2x2), for w = [w_r; w_b].
%     Rd       - Measurement noise covariance (1x1), variance of psi_meas.
%     psi_meas - Noisy yaw-angle measurement at step k [rad].
%     delta    - Rudder angle applied at step k [rad] (filter input).
%
%   OUTPUTS
%     x_pst    - Corrected (a posteriori) estimate at step k, used as
%                feedback in the heading autopilot.
%     P_pst    - Corrected error covariance at step k.
%     x_prd    - Predicted state for step k+1.
%     P_prd    - Predicted error covariance for step k+1.
%     innov    - Innovation at step k [rad]: wrapped measurement residual
%                y - Cd x_prd, computed from the PREDICTED state.
%     S        - Innovation covariance at step k (1x1). Together with
%                innov it gives the consistency test NIS = innov^2 / S.

    I = eye(size(P_prd));

    % Corrector: innovation, Kalman gain and measurement update
    innov = ssa(psi_meas - Cd*x_prd);              % wrapped innovation [rad]
    S     = Cd*P_prd*Cd' + Rd;                     % innovation covariance
    K_kf  = P_prd*Cd' / S;                         % Kalman gain
    x_pst = x_prd + K_kf*innov;
    x_pst(1) = ssa(x_pst(1));                      % keep yaw estimate in (-pi, pi]
    P_pst = (I - K_kf*Cd)*P_prd*(I - K_kf*Cd)' + K_kf*Rd*K_kf';   % Joseph form

    % Predictor: propagate to the next time step
    x_prd = Ad*x_pst + Bd*delta;
    P_prd = Ad*P_pst*Ad' + Ed*Qd*Ed';
end
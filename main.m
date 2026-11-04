% Project in TTK4190 Guidance, Navigation and Control of Vehicles 
%
% Author:           Antonio Zerbato
% Study program:    Robotics Engineering 

% Add folder for 3-D visualization files
addpath(genpath('flypath3d_v3'))
addpath(genpath('pathplotter and waypoints'))
clear WP_selector      % reset persistent segment index
clear ILOS_guidance    % reset persistent segment index

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
T_final = 6500;	        % Final simulation time (s)
h = 0.1;                % Sampling time (s)

psi_ref = 10 * pi/180;  % desired yaw angle (rad)
U_ref   = 9;            % desired surge speed (m/s) 
Vc = 1;                 %current vel 
beta_Vc = deg2rad(45);  %current angle

Vw = 10;                %wind vel
beta_Vw = deg2rad(135);
rho_a   = 1.247;
cy      = 0.95;
cn      = 0.15;
L_oa    = 161;          %ship lenght      
A_Lw    = 10 * L_oa;
t_wind  = 200;          % wind starts after 200 s
delta_max  = 40 * pi/180;        % max rudder angle      (rad)
current_on  = true;
crab_compensation = true;
feedback = 'kf';



% initial states
eta_0 = [0 0 -110*pi/180]';
nu_0  = [0 0 0]';
delta_0 = 0;
n_0 = 0;
Qm_0    = 0;
x = [nu_0' eta_0' delta_0 n_0 Qm_0]'; % The state vector can be extended with addional states here

%% Nomoto Model
u_design = 9;     % speed used to linearize the model for the heading autopilot
[T,K] = nomoto_first_order(u_design);

%% gain for PID
wb = 0.005;      %control bandwidth (inital 0.06)
zeta = 1;       %damping ratio
[Kp, Kd, Ki] = pid_desing_nomoto(wb, zeta,T, K );
Ki= 0;
xd = [eta_0(3); 0; 0];          % reference model states: start at initial heading, at rest
e_int = 0;                       % PID integral state


%% State space system
[Ad, Bd, Cd, Ed] = discrete_system(K,T, h);
fprintf('rank(obsv) = %d\n', rank(obsv(Ad,Cd))); 

% Sensor noise (Problem 4a)
sigma_psi = deg2rad(0.5);                 % yaw angle std [rad]
sigma_r   = deg2rad(0.1);                 % yaw rate std  [rad/s]

% Covariances (tune these, see section 4)
Rd    = sigma_psi^2;                      % measurement noise (yaw only)
Qd    = diag([1e-6, 1e-8]);               % process noise [w_r, w_b]
P_prd = diag([pi^2, deg2rad(1)^2, deg2rad(5)^2]);   % initial covariance P(0)
x_prd = [0; 0; 0];                        % initial estimate [psi; r; b]




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:h:T_final;                % Time vector
nTimeSteps = length(t);         % Number of time steps

simdata = zeros(nTimeSteps, 24); % Pre-allocate matrix for efficiency

for i = 1:nTimeSteps

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 3, 4a) Noisy measurements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    psi_meas = ssa(x(6) + sigma_psi*randn);   % gyrocompass
    r_meas   = x(3) + sigma_r*randn;          % rate gyro

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 3, 4b) Kalman filter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x_pst,P_pst,x_prd,P_prd,innov,S] = KF(x_prd,P_prd,Ad,Bd,Ed,Cd,Qd,Rd,psi_meas,x(7));
    NIS   = innov^2 / S;                              % ~ chi2 with 1 dof
    e_est = [ssa(x_pst(1) - x(6)); x_pst(2) - x(3)];  % estimation error (psi, r)
    NEES  = e_est' / P_pst(1:2,1:2) * e_est;          % ~ chi2 with 2 dof


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 3, 4c/4d) Feedback signals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch feedback
        case 'true',  psi_fb = x(6);      r_fb = x(3);
        case 'noisy', psi_fb = psi_meas;  r_fb = r_meas;
        case 'kf',    psi_fb = x_pst(1);  r_fb = x_pst(2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 1a) Add current here 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if current_on
        uc = Vc * cos(beta_Vc - x(6));
        vc = Vc * sin(beta_Vc - x(6));
    else 
        uc = 0;
        vc = 0;
    end
    nu_c = [ uc vc 0 ]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 3, 1a) Guidance Law
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [chi_d, last, e_y] = guidance(x(4), x(5), h);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 3, 2c) Crab compensation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if crab_compensation
        U = sqrt(x(1)^2 + x(2)^2);              % speed over ground
        if U > 1                                % prevents the intial spike
            beta_c = atan2(x(2), x(1));         % crab angle (velocity over ground)
        else
            beta_c = 0;
        end
        psi_ref = ssa(chi_d - beta_c);
    else
        psi_ref = chi_d;
    end
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 1c) Add wind here 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if t(i) > t_wind
        u_w  = Vw * cos(beta_Vw - x(6));
        v_w  = Vw * sin(beta_Vw - x(6));
        u_rw = x(1) - u_w;
        v_rw = x(2) - v_w;
        V_rw = sqrt(u_rw^2 + v_rw^2);
        gamma_rw = -atan2(v_rw, u_rw);
        q     = 0.5 * rho_a * V_rw^2;
        Ywind = q * cy * sin(gamma_rw)   * A_Lw;
        Nwind = q * cn * sin(2*gamma_rw) * A_Lw * L_oa;
    else
        Ywind = 0;
        Nwind = 0;
    end
    tau_wind = [0 Ywind Nwind]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 2d) Add a reference model here 
    % Define it as a function
    % check eq. (15.143) in (Fossen, 2021) for help
    %
    % The result should look like this:
    % xd_dot = ref_model(xd,psi_ref(i));
    % psi_d = xd(1);
    % r_d = xd(2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if t(i) > 400 %comment to compute point 3e (constant heading) 
    %     psi_ref = -20* pi/180;
    % end
    xd_dot = ref_model(xd,psi_ref );
    psi_d = xd(1);
    r_d = xd(2);
    u_d = U_ref;



    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 2d) Add the heading controller here 
    % Define it as a function
    %
    % The result should look like this:
    % delta_c =  PID_heading(e_psi,e_r,e_int);
    e_psi = ssa(psi_fb - psi_d);                              %psi error
    e_r = r_fb - r_d;                                       %r error
    delta_c = PID_heading(e_psi,e_r,e_int, Kp,Kd, Ki);      %rudder angle command

    saturation = abs(delta_c) > delta_max;
    if ~saturation
     e_int = e_int + h*e_psi;                            %update integral term
    end 
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delta_c = 0.1;              % rudder angle command (rad)
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 3e) Add open loop speed control here
    % Define it as a function
    %
    % The result should look like this:
    % n_c = open_loop_speed_control(U_ref);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % n_c = 10;                   % propeller speed [radians per second (rps)]
    n_c = open_loop_speed_control(U_ref);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 3f) Replace the open loop speed controller, 
    % with a closed loop speed controller here 
    % Define it as a function
    %
    % The result should look like this:
    % n_c = closed_loop_speed_control(u_d,e_u,e_int_u);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ship dynamics
    u = [delta_c n_c]';
    [xdot,tau_total] = ship(x,u,nu_c,tau_wind);
    
    % store simulation data in a table (for testing)
    simdata(i,:) = [x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) u_d psi_d r_d uc vc e_y chi_d ...
                    psi_meas r_meas x_pst' NIS NEES];
   
    % Euler integration
    % x = euler2(xdot,x,h); 
    % Runge Kutta 4 integration
    x = rk4(@ship,h,x,u,nu_c,tau_wind);
    xd = xd + h*xd_dot;                     %integration of reference model

    % --- Progress Update Logic ---
    if mod(i, floor(nTimeSteps / 10)) == 0 % Print an update every 10%
        progress = (i/nTimeSteps) * 100;
        fprintf('  %d%% complete\n', round(progress));
    elseif i == 1
        disp("  Simulation in progress:")
        fprintf('  %d%% complete\n', 0);
    end

end
simdata = simdata(1:i,:);
t = t(1:i);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u           = simdata(:,1);                 % m/s
v           = simdata(:,2);                 % m/s
r           = simdata(:,3);                 % rad/s
r_deg       = (180/pi) * r;                 % deg/s
x           = simdata(:,4);                 % m
y           = simdata(:,5);                 % m
psi         = simdata(:,6);                 % rad
psi_deg     = (180/pi) * psi;               % deg
delta_deg   = (180/pi) * simdata(:,7);      % deg
n           = (30/pi) * simdata(:,8);       % rpm
delta_c_deg = (180/pi) * simdata(:,9);      % deg
n_c         = (30/pi) * simdata(:,10);      % rpm
u_d         = simdata(:,11);                % m/s
psi_d       = simdata(:,12);                % rad
psi_d_deg   = (180/pi) * psi_d;             % deg
r_d         =  simdata(:,13);               % rad/s
r_d_deg     = (180/pi) * r_d;               % deg/s
uc = simdata(:, 14);
vc = simdata(:, 15);
e_y = simdata(:, 16);
chi_d = (180/pi)* simdata(:, 17);
psi_meas = simdata(:,18);   r_meas = simdata(:,19);
psi_hat  = simdata(:,20);   r_hat  = simdata(:,21);   b_hat = simdata(:,22);
NIS_log  = simdata(:,23);
NEES_log = simdata(:,24);
u_r = u - uc;
v_r = v - vc;

%% PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter consistency: NIS (1 dof) and NEES (2 dof)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx   = t > 100;                 % skip initial convergence
tt    = t(idx);
nis   = NIS_log(idx);
nees  = NEES_log(idx);
alpha = [0.025 0.975];           % 95% two-sided interval

% Per-sample 95% bounds of the chi-square distribution
if exist('chi2inv','file')
    B_NIS  = chi2inv(alpha, 1);
    B_NEES = chi2inv(alpha, 2);
else                             % no Statistics Toolbox
    B_NIS  = [0.00098 5.024];
    B_NEES = [0.0506  7.378];
end

ANIS    = mean(nis);
ANEES   = mean(nees);
in_NIS  = 100*mean(nis  >= B_NIS(1)  & nis  <= B_NIS(2));
in_NEES = 100*mean(nees >= B_NEES(1) & nees <= B_NEES(2));
fprintf('ANIS  = %.3f (ideal 1), %.1f%% of samples inside 95%% bounds\n', ANIS,  in_NIS);
fprintf('ANEES = %.3f (ideal 2), %.1f%% of samples inside 95%% bounds\n', ANEES, in_NEES);
figure(10); clf

% ---- NIS ----
subplot(211); hold on; grid on
plot(tt, nis, '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 3);
yline(B_NIS(2), 'r--', 'linewidth', 1.5);
yline(B_NIS(1), 'r--', 'linewidth', 1.5);
yline(ANIS,     'b-',  'linewidth', 2);
set(gca, 'YScale', 'log'); ylim([1e-4 50])
title(sprintf('NIS: average %.2f (ideal 1), %.1f%% inside 95%% bounds', ANIS, in_NIS));
xlabel('Time (s)'); ylabel('NIS');
legend('NIS', '97.5% bound', '2.5% bound', 'average', 'Location', 'eastoutside')

% ---- NEES ----
subplot(212); hold on; grid on
plot(tt, nees, '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 3);
yline(B_NEES(2), 'r--', 'linewidth', 1.5);
yline(B_NEES(1), 'r--', 'linewidth', 1.5);
yline(ANEES,     'b-',  'linewidth', 2);
set(gca, 'YScale', 'log'); ylim([1e-3 100])
title(sprintf('NEES: average %.2f (ideal 2), %.1f%% inside 95%% bounds', ANEES, in_NEES));
xlabel('Time (s)'); ylabel('NEES');
legend('NEES', '97.5% bound', '2.5% bound', 'average', 'Location', 'eastoutside')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 4) State estimation plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi_true_deg = rad2deg(ssa(psi));          % wrap true yaw for comparison
psi_meas_deg = rad2deg(psi_meas);
psi_hat_deg  = rad2deg(psi_hat);
r_meas_deg   = rad2deg(r_meas);
r_hat_deg    = rad2deg(r_hat);

% 4a) Noisy measurements vs true states
figure(6); clf
subplot(211)
plot(t, psi_meas_deg, t, psi_true_deg, 'linewidth', 1.5); grid on
title('Yaw angle: measurement vs true'); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('measured \psi', 'true \psi')
subplot(212)
plot(t, r_meas_deg, t, r_deg, 'linewidth', 1.5); grid on
title('Yaw rate: measurement vs true'); xlabel('Time (s)'); ylabel('Rate (deg/s)');
legend('measured r', 'true r')

% 4b) Kalman filter estimates vs true states, and bias estimate
figure(7); clf
subplot(311)
plot(t, psi_true_deg, t, psi_hat_deg, '--', 'linewidth', 1.5); grid on
title('Yaw angle: estimate vs true'); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('true \psi', 'estimated \psi')
subplot(312)
plot(t, r_deg, t, r_hat_deg, '--', 'linewidth', 1.5); grid on
title('Yaw rate: estimate vs true'); xlabel('Time (s)'); ylabel('Rate (deg/s)');
legend('true r', 'estimated r')
subplot(313)
plot(t, rad2deg(b_hat), 'linewidth', 1.5); grid on
title('Estimated rudder bias'); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('estimated b')

% 4b) Estimation errors (zoomed, after initial convergence)
idx = t > 50;
figure(8); clf
subplot(211)
plot(t(idx), rad2deg(ssa(psi_hat(idx) - psi(idx))), 'linewidth', 1.2); grid on
title('Yaw angle estimation error'); xlabel('Time (s)'); ylabel('\psi_{hat} - \psi (deg)');
subplot(212)
plot(t(idx), r_hat_deg(idx) - r_deg(idx), 'linewidth', 1.2); grid on
title('Yaw rate estimation error'); xlabel('Time (s)'); ylabel('r_{hat} - r (deg/s)');

% 4c/4d) Closed-loop performance with the selected feedback
figure(9); clf
sgtitle(['Feedback: ' feedback])
subplot(221)
plot(y, x, 'linewidth', 1.5); axis equal; grid on
title('North-East position'); xlabel('East (m)'); ylabel('North (m)');
subplot(222)
plot(t, delta_c_deg, 'linewidth', 1.2); grid on
ylim([-1.2 1.2]*rad2deg(delta_max))
title('Desired rudder angle'); xlabel('Time (s)'); ylabel('\delta_c (deg)');
subplot(223)
plot(t, psi_true_deg, t, psi_d_deg, '--', 'linewidth', 1.5); grid on
title('Yaw angle'); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('\psi', '\psi_d')
subplot(224)
plot(t, r_deg, t, r_d_deg, '--', 'linewidth', 1.5); grid on
title('Yaw rate'); xlabel('Time (s)'); ylabel('Rate (deg/s)');
legend('r', 'r_d')



        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3, 2a) Side Slip and Crab Angle
beta_crab  = atan2(v, u);      % crab angle (over ground)
beta_slip  = atan2(v_r, u_r);    % sideslip angle (through water)

chi = rad2deg(ssa(psi + beta_crab));
figure(4)
figure(gcf)
subplot(211)
plot(t, rad2deg(beta_slip), t, rad2deg(beta_crab), 'linewidth', 2);
title('Sideslip and crab angle'); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('sideslip \beta', 'crab \beta_c')
subplot(212)
plot(t,e_y, 'linewidth', 2);
title('Cross track error'); xlabel('Time (s)'); ylabel('Error (m)');
legend('cross track error')


figure(5)
figure(gcf)
subplot(212)
plot(t, chi_d, t, chi, 'linewidth', 2);
title('Course desired and course'); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('course desired \chi_d', 'course \chi')
subplot(211)
plot(t, chi_d - chi,'linewidth', 2);
title("Course error"); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('course error \e')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





figure(3)
figure(gcf)
subplot(311)
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions'); xlabel('(m)'); ylabel('(m)'); 
subplot(312)
plot(t,psi_deg,t,psi_d_deg,'linewidth',2);
title('Actual and desired yaw angle'); xlabel('Time (s)');  ylabel('Angle (deg)'); 
legend('actual yaw','desired yaw')
subplot(313)
plot(t,r_deg,t,r_d_deg,'linewidth',2);
title('Actual and desired yaw rates'); xlabel('Time (s)');  ylabel('Angle rate (deg/s)'); 
legend('actual yaw rate','desired yaw rate')

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocity'); xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('actual surge','desired surge')
subplot(312)
plot(t,n,t,n_c,'linewidth',2);
title('Actual and commanded propeller speed'); xlabel('Time (s)'); ylabel('Motor speed (RPM)');
legend('actual RPM','commanded RPM')
subplot(313)
plot(t,delta_deg,t,delta_c_deg,'linewidth',2);
title('Actual and commanded rudder angle'); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('actual rudder angle','commanded rudder angle')
%% Create objects for 3-D visualization 
% Since we only simulate 3-DOF we need to construct zero arrays for the 
% excluded dimensions, including height, roll and pitch
z = zeros(length(x),1);
phi = zeros(length(psi),1);
theta = zeros(length(psi),1);

% create object 1: ship (ship1.mat)
new_object('flypath3d_v3/ship1.mat',[x,y,z,phi,theta,psi],...
'model','royalNavy2.mat','scale',(max(max(abs(x)),max(abs(y)))/1000),...
'edge',[0 0 0],'face',[0 0 0],'alpha',1,...
'path','on','pathcolor',[.89 .0 .27],'pathwidth',2);
% 
% % Plot trajectories 
% flypath('flypath3d_v3/ship1.mat',...
% 'animate','on','step',500,...
% 'axis','on','axiscolor',[0 0 0],'color',[1 1 1],...
% 'font','Georgia','fontsize',12,...
% 'view',[-25 35],'window',[900 900],...
% 'xlim', [min(y)-0.1*max(abs(y)),max(y)+0.1*max(abs(y))],... 
% 'ylim', [min(x)-0.1*max(abs(x)),max(x)+0.1*max(abs(x))], ...
% 'zlim', [-max(max(abs(x)),max(abs(y)))/100,max(max(abs(x)),max(abs(y)))/20]); 


%Plotter
pathplotter(x,y);


V_rated = 750; % V
Ra = 0.044;% Ohm
La = 0.000748; % H
Ke = 0.3255; % V*s /rad
Kt = Ke; % torque constant
T_rated = 277.5; % Nm
J = 0.230000; % kg*m^2
n_rated = 19800; % rev/ min

%% calculate the controller parameters
Ta = La/Ra; % electrical time constant
Tm = J*Ra/ (Kt*Ke); % mechanical time constant
Kob = V_rated/Ke; % static gain
Kp = Tm/(2*Ta*Kob); % proportional coefficient
Ti = Tm;% integral time constant
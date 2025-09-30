% Transformer nameplate information: 37.5 kVA, 7200/240 V, 60 Hz
% Add following parameters to m-file: rated primary voltage, rated secondary voltage,
% rated primary current, rated secondary current.
rpv=7200;
rsv=240;
rpc=37500/7200;
rsc=37500/240;
%% Open-circuit test
% Simulate the open-circuit test and then add following parameters to m-file: Voc, Ioc, Poc
V1oc = 7200;
V2oc = 242.7;
Ioc = 0.6489;
Poc = 600.8;

% Write script for determination of the magnetization branch parameters (Rc, Xm)
a = V1oc/V2oc; %a=30
Zoc = V1oc/Ioc; %11096
PF = Poc/(V1oc*Ioc); %0.1286
Rc = Zoc * PF; %1426.8
Xm = Zoc * sqrt(1-PF^2); %11004
theta_oc = acos(Poc/(V1oc*Ioc)); %1.4418

%draw
vVocdZ = V1oc / sqrt(Xm^2+Rc^2) * exp(1i * 0);
vIoc = Ioc * exp(1i * theta_oc);
vIocL = V1oc / Xm * exp(1i * -pi/2);
vIocR = V1oc / Rc * exp(1i * 0);
figure;
hold on;
quiver(0, 0, real(vVocdZ), imag(vVocdZ), 0, 'cyan', 'LineWidth', 3);
quiver(0, 0, real(vIoc), imag(vIoc), 0, 'green', 'LineWidth', 2);
quiver(0, 0, real(vIocL), imag(vIocL), 0, 'blue', 'LineWidth', 2);
quiver(0, 0, real(vIocR), imag(vIocR), 0, 'magenta', 'LineWidth', 2);
hold off;
axis equal;
title("Open Circuit Test Phasor Diagram");
legend("Voc/Z", "Ioc", "IocL", "IocR");
saveas(gcf, "Open_Circuit_Test_Phasor_Diagram.png");


%% Short - circuit test
% Simulate the short-circuit test and then add following parameters to m-file: Vsc, Isc, Psc
Vsc = 360;
Isc = 7.045;
Psc = 694.4;

% Write script for determination of the leakage resistances (R1, R2) and reactances (Xl1, Xl2) 
PF_sc = Psc/(Vsc*Isc);
Z_eq = Vsc/Isc;
R_eq = Psc/Isc^2;
X_eq = sqrt(Z_eq^2-R_eq^2);
a = 7200/240;
R1 = R_eq/2;%6.9955
R2 = R_eq/(2*a^2)%0.0077
X11 = X_eq/2;%24.5737
X12 = X_eq/(2*a^2)%0.0273

%draw
theta_sc = acos(PF_sc);

vVsc = Vsc * exp(1i*0);
vVscR = Isc * R_eq *exp(-1i * theta_sc);
vVscX = Isc * X_eq * exp(1i * (pi/2-theta_sc));
figure;
hold on;
quiver(0, 0, real(vVsc), imag(vVsc), 0, 'r', 'LineWidth', 2);
quiver(0, 0, real(vVscR), imag(vVscR), 0, 'green', 'LineWidth', 2);
quiver(real(vVscR), imag(vVscR),real(vVscX),imag(vVscX),0, 'blue', 'LineWidth', 2);
axis equal;
title("Short Circuit Test Phasor Diagram");
legend("Vsc", "IscRsc", "IscXsc");
saveas(gcf,"Short_Circuit_Test_Phasor_Diagram.png");
%% Load tests
% R-Load tests
% Values for R-Load:
R1 = [1.4 1.5 1.6 1.7 2 3 4 5 7 10]; % ohm
% Simulate the model with each value of the R-load and then add the
% following parameters to m-file for each test:
U1_R = [7200 7200 7200 7200 7200 7200 7200 7200 7200 7200];
U2_R = [239.8 240 240.1 240.3 240.7 241.4 241.7 241.9 242.1 242.3];
I1_R = [5.915 5.536 5.204 4.91 4.206 2.879 2.223 1.837 1.409 1.109];
I2_R = [171.3 160 150.1 141.4 120.3 80.45 60.42 48.38 34.59 24.23];
P1_R = [42130 39390 37000 3488029790 20120 15260 12340 8994 6480];
P2_R = [41060 38380 36040 33960 28960 19420 14600 11700 8374 5870];
% Write script for determination of the load ratio, transformer efficiency, power
beta = zeros( size(R1));
power_factor = zeros(size(R1));
efficiency = zeros(size(R1));
voltage_regulation_mod = zeros(size(R1));
voltage_regulation_calc = zeros(size(R1));
theta_r = zeros(size(R1));
for k = 1:length(R1)
    beta(k) = I2_R(k)/rsc ; 
    power_factor(k) = P2_R(k)/(U2_R(k)*I2_R(k)); 
    efficiency(k) = P2_R(k) /(P2_R(k) + beta(k)^2 *R_eq *I1_R(k)^2 + Poc); 
    theta_r(k) = acos(power_factor(k)); % Calculate phase angle
    voltage_regulation_mod(k) = ((V2oc - U2_R(k))* 100 )/ U2_R(k); % Calculate measured voltage regulation rate
    voltage_regulation_calc(k)= (beta(k) *(Vsc/a) * (cos(theta_sc) * cos(theta_r(k)) + sin(theta_sc) * sin(theta_r(k)))) *100 /U2_R(k);
end
    
% factor, voltage regulation

%draw
vV2R = U2_R(1)* exp(1i* 0);
S_1_r = U1_R(1)*I1_R(1);
PF_1_r = P1_R(1)/S_1_r;
theta_prime_r = acos(PF_1_r);

vV1R = U1_R(1)* exp(1i* theta_prime_r);
vV2to1R = vV2R * a;

vReqI1R = I1_R(1)* (R_eq + R1(1)/ a^2)* exp(1i* 0)*exp(-1i * theta_r(1));
vXeqI1R = I1_R(1)* X_eq * exp(1i * pi/2)*exp(-1i * theta_r(1));
quiver(0,0, real(vV1R), imag(vV1R),0, 'r', 'Linewidth',2);
hold on;
quiver(0,0,real(vV2to1R),imag(vV2to1R),0,'g', 'Linewidth', 2);
quiver(real(vV2to1R),imag(vV2to1R),real(vReqI1R),imag(vReqI1R),0,'b', 'Linewidth', 2);
quiver(real(vReqI1R + vV2to1R),imag(vReqI1R + vV2to1R),real(vXeqI1R),imag(vXeqI1R),0,'m', 'Linewidth', 2);
axis equal;
title("R Load Test Phasor Diagram");
legend("v1","V2to1","ReqI1","XeqI1");

% Draw diagrams
figure;
subplot(3, 2, 1);
plot(beta, P2_R, 'r');
xlabel('Load Ratio (β)');
ylabel('P2 (W)');
title('P2 vs. β');
grid on;

subplot(3, 2, 2);
plot(beta, power_factor, 'r');
xlabel('Load Ratio (β)');
ylabel('Power Factor (PF_{load})');
title('Power Factor vs. β');
grid on;

subplot(3, 2, 3);
plot(beta, efficiency, 'r');
xlabel('Load Ratio (β)');
ylabel('Efficiency (η)');
title('Efficiency vs. β');
grid on;

figure;
plot(beta, voltage_regulation_calc, 'r', 'DisplayName', 'Calculated');
hold on;
plot(beta, voltage_regulation_mod, 'b', 'DisplayName', 'Measured');
xlabel('Load Ratio (β)');
ylabel('Voltage Regulation (ΔV%)');
title('Voltage Regulation vs. β');
legend('show');
grid on;
hold off;





sgtitle('Performance Characteristics');



% RL-Load tests
% Values for R-Load:
R2 = [1.23 1.32 1.41 1.49 1.76 2.64 3.52 4.4 6.16 8.8]; % ohm
L = [1.76 1.89 2.01 2.14 2.52 3.78 5.04 6.3 8.82 12.6]; %mH
% Simulate the model with each value of the RL-load and then add the
% following parameters to m-file for each test:
U1_RL=[7200 7200 7200 7200 7200 7200 7200 7200 7200 7200];
U2_RL=[235.7 236.2 236.6 236.9 237.8 239.4 240.2 240.7 241.3 241.7];
I1_RL=[6.1 5.723 5.398 5.131 4.427 3.119 2.464 2.073 1.63 1.306];
I2_RL=[168.7 157.5 147.8 139.8 118.9 79.8 60.05 48.14 34.47 24.17];
P1_RL=[36060 33730 31760 30050 25710 17520 13360 10840 7940 5752];
P2_RL=[34990 32720 30800 29130 24880 16810 12690 10200 7317 5140];
% Write script for determination of the load ratio, transformer efficiency, power
beta_RL = zeros( size(R2));
PF_RL = zeros(size(R2));
theta_RL = zeros(size(R2));
voltage_regulation_mod_RL = zeros(size(R2));
voltage_regulation_calc_RL = zeros(size(R2));
eta_RL = zeros(size(R2));
for k = 1:length(R2)
    beta_RL(k) = I2_RL(k)/rsc ; 
    PF_RL(k) = P2_RL(k)/(U2_RL(k)*I2_RL(k)); 
    eta_RL(k) = P2_RL(k) /(P2_R(k) + I1_RL(k)^2 *R_eq *beta_RL(k)^2 + Poc); 
    theta_RL(k) = acos(PF_RL(k)); % Calculate phase angle
    voltage_regulation_mod_RL(k) = ((V2oc - U2_RL(k))* 100 )/ U2_RL(k); % Calculate measured voltage regulation rate
    voltage_regulation_calc_RL(k)= (beta_RL(k) *(Vsc/a) * (cos(theta_sc) * cos(theta_RL(k)) + sin(theta_sc) * sin(theta_RL(k)))) *100 /U2_RL(k);
end
% factor, voltage regulation

%draw
vV2RL = U2_RL(1)* exp(1i* 0);
S_1_rL = U1_RL(1)*I1_RL(1);
PF_1_rL = P1_RL(1)/S_1_rL;
theta_prime_rL = acos(PF_1_rL);

vV1RL = U1_RL(1)* exp(1i* theta_prime_rL- theta_RL(1));
vV2to1RL = vV2RL * a;

vReqI1RL = I1_RL(1)* (R_eq + R2(1)/ a^2)* exp(1i* 0)*exp(-1i * theta_RL(1));
vXeqI1RL = I1_RL(1)* X_eq * exp(1i * pi/2)+(60*2*pi*L(1)/1000/a^2*exp(1i*pi/2))*exp(-1i * theta_r(1));
quiver(0,0, real(vV1RL), imag(vV1RL),0, 'r', 'Linewidth',2);
hold on;
quiver(0,0,real(vV2to1RL),imag(vV2to1RL),0,'g', 'Linewidth', 2);
quiver(real(vV2to1RL),imag(vV2to1RL),real(vReqI1RL),imag(vReqI1RL),0,'b', 'Linewidth', 2);
quiver(real(vReqI1RL + vV2to1RL),imag(vReqI1RL + vV2to1RL),real(vXeqI1RL),imag(vXeqI1RL),0,'m', 'Linewidth', 2);
hold off;
axis equal;
title("RL Load Test Phasor Diagram");
legend("v1","V2to1","ReqI1","XeqI1");


figure;
subplot(3, 2, 1);
plot(beta_RL, P2_RL, 'r');
xlabel('β');
ylabel('P2 (W)');
title('P2 vs β');
grid on;

subplot(3, 2, 2);
plot(beta_RL, PF_RL, 'r');
xlabel('β');
ylabel('PFload');
title('PFload vs β');
grid on;

subplot(3, 2, 3);
plot(beta_RL, eta_RL, 'r');
xlabel('β');
ylabel('η');
title('η vs β');
grid on;

figure;
plot(beta_RL, voltage_regulation_calc_RL, 'r', 'DisplayName', 'ΔV%2_calc');
hold on;
plot(beta_RL, voltage_regulation_mod_RL, 'b', 'DisplayName', 'ΔV%2_mod');
xlabel('β');
ylabel('ΔV%');
title('ΔV%2_calc vs ΔV%2_mod vs β');
legend('show');
grid on;
hold off;




% RC-Load tests
% Values for RC-Load:
R3 = [1.05 1.125 1.2 1.27 1.5 2.25 3 3.75 5.25 7.5]; %ohm
C = [2.86 2.67 2.51 2.36 2.01 1.34 1 0.8 0.57 0.4]; %mF
% Simulate the model with each value of the RC-load and then add the
% following parameters to m-file for each test:
U1_RC=[7200 7200 7200 7200 7200 7200 7200 7200 7200 7200];
U2_RC = [246.9 246.7 246.4 246.2 245.7 244.7 244.2 243.9 243.6 243.3];
I1_RC= [5.626 5.224 4.881 4.58 3.836 2.46 1.783 1.394 0.9731 0.7076];
I2_RC = [176.3 164.3 154.1 145.2 123 81.65 60.98 48.73 34.72 24.3];
P1_RC= [33700 31400 29460 27690 23510 15700 11810 9535 6943 5036];
P2_RC =[32620 30380 28490 26770 22680 15000 11150 8902 6327 4429];
% Write script for determination of the load ratio, transformer efficiency, power
beta_RC = zeros( size(R3));
PF_RC = zeros(size(R3));
theta_RC = zeros(size(R3));
voltage_regulation_mod_RC = zeros(size(R3));
voltage_regulation_calc_RC = zeros(size(R3));
eta_RC = zeros(size(R3));
for k = 1:length(R3)
    beta_RC(k) = I2_RC(k)/rsc ; 
    PF_RC(k) = P2_RC(k)/(U2_RC(k)*I2_RC(k)); 
    eta_RC(k) = P2_RC(k) /(P2_R(k) + I1_RC(k)^2 *R_eq *beta_RC(k)^2 + Poc); 
    theta_RC(k) = acos(PF_RC(k)); % Calculate phase angle
    voltage_regulation_mod_RC(k) = ((V2oc - U2_RC(k))* 100 )/ U2_RC(k); % Calculate measured voltage regulation rate
    voltage_regulation_calc_RC(k)= (beta_RC(k) *(Vsc/a) * (cos(theta_sc) * cos(theta_RC(k)) - sin(theta_sc) * sin(theta_RC(k)))) *100 /U2_RC(k);
    
end


figure;
subplot(3, 2, 1);
plot(beta_RC, P2_RC, 'r');
xlabel('β');
ylabel('P2 (W)');
title('P2 vs β');
grid on;

subplot(3, 2, 2);
plot(beta_RC, PF_RC, 'r');
xlabel('β');
ylabel('PFload');
title('PFload vs β');
grid on;

subplot(3, 2, 3);
plot(beta_RC, eta_RC, 'r');
xlabel('β');
ylabel('η');
title('η vs β');
grid on;



% Plotting voltage_regulation_calc_RC and voltage_regulation_mod_RC on the same graph
figure;
plot(beta_RC, voltage_regulation_calc_RC, 'r', 'DisplayName', 'ΔV%2_calc');
hold on;
plot(beta_RC, voltage_regulation_mod_RC, 'b', 'DisplayName', 'ΔV%2_mod');
xlabel('β');
ylabel('ΔV%');
title('ΔV%2_calc vs ΔV%2_mod vs β');
legend('show');
grid on;
hold off;
% factor, voltage regulation
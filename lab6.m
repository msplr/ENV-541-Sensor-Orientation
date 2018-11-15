%% data loading
[data, fIMU] = readimu('1109_1522_PostProBinaryDecoded.imu', 'IXSEA');
t_start = 485090;
t_stop = 485200;

filter = data(:,1) > t_start & data(:,1) < t_stop;
t = data(filter,1); t = t - t(1);
gyro = data(filter,2:4);
accel = data(filter,5:7);

% convert NWU to NED frame
R_NNW2NED = [1 0 0; 0 -1 0; 0 0 -1];
gyro = (R_NNW2NED*gyro')';
accel = (R_NNW2NED*accel')';

%% Reference comparison
% reference values
latitude_EPFL = 46+(31+17/60)/60; % North [deg]
w_ref = 7.2921150e-5; % [rad/s]
g_EPFL = 9.8055; % [m/s^2]

% Reference from IXSEA software in NWD
latitude_ref = 46.52094807; % [deg]
heading_ref = 41.978; % [deg]
roll_ref = -2.653; % [deg]
pitch_ref = -2.464; % [deg]

% convert reference NWD -> NED
pitch_ref = -pitch_ref;

accel = mean(accel,1)';
gyro = mean(gyro,1)';
gyro_norm = vecnorm(gyro);
accel_norm = vecnorm(accel);

fprintf('Sensor norm:\n');
fprintf('gyro: %.7e rad/s, difference to w_ref %.3e rad/s\n', gyro_norm, gyro_norm-w_ref);
fprintf('accel: %.7e m/s^2, difference to g_EPFL %.3e m/s^2\n', accel_norm, accel_norm-g_EPFL);

%% Test (simulate sensor readings)
if 0
    r = 23 *pi/180;
    p = 42 *pi/180;
    phi = 56 *pi/180;
    y = 78 *pi/180;
    RR = R1(r)*R2(p);
    accel = RR*[0;0;-g_EPFL];
    gyro = RR*R3(y)*[w_ref*cos(phi);0 ;-w_ref*sin(phi)];
    gyro_norm = vecnorm(gyro);
    accel_norm = vecnorm(accel);
end

%% Accelerometer leveling
roll = atan(accel(2)/accel(3));
pitch = asin(accel(1)/accel_norm);
fprintf('roll %.3f deg, pitch %.3f deg\n', roll*180/pi, pitch*180/pi);
fprintf('reference: roll %.3f deg, pitch %.3f deg\n', roll_ref, pitch_ref);

%% Gyrocompassing
Rlb = R1(roll)*R2(pitch);
Rbl = Rlb';
wl = Rbl * gyro;

yaw = atan(-wl(2)/wl(1));
fprintf('yaw %.3f deg\n', yaw*180/pi);
fprintf('reference yaw %.3f deg\n', heading_ref);

latitude = asin(-wl(3)/gyro_norm);
fprintf('latitude %.3f deg\n', latitude*180/pi);
fprintf('reference latitude %.3f deg\n', latitude_ref);

%% functions
function [R] = R1(a)
    R = [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)];
end

function [R] = R2(a)
    R = [cos(a) 0 -sin(a); 0 1 0; sin(a) 0 cos(a)];
end

function [R] = R3(a)
    R = [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];
end

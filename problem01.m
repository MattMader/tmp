%% Setup

% clear
clear
clc

% Load functions
addpath functions

%% Part a

% Euler (3-2-1) attitude at t=0
theta1 = -pi/4; % [rad]
theta2 = pi/8; % [rad]
theta3 = pi/5; % [rad]

% DCM attitude at t=0
C = EA321toDCM(theta1,theta2,theta3); % [-]

% Euler Parameter attitude at t=0
epsilon0 = DCMtoEP(C); % [-]

% Confirm unit-norm
disp(norm(epsilon0))

%% Part c

% time span
tspan = 0:0.5:100; % [s]

% integration options
opts = odeset("RelTol",1e-10,"AbsTol",1e-10);

% integration
[t,epsilon] = ode45(@odefun1,tspan,epsilon0,opts);

% plots
figure(1)
subplot(2,2,1)
plot(t,epsilon(:,1),'b')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Euler Parameter $\epsilon_1$","Interpreter","latex")
title("$\epsilon_1(t)$","Interpreter","latex")

subplot(2,2,2)
plot(t,epsilon(:,2),'b')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Euler Parameter $\epsilon_2$","Interpreter","latex")
title("$\epsilon_2(t)$","Interpreter","latex")

subplot(2,2,3)
plot(t,epsilon(:,3),'b')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Euler Parameter $\epsilon_3$","Interpreter","latex")
title("$\epsilon_3$(t)","Interpreter","latex")

subplot(2,2,4)
plot(t,epsilon(:,4),'b')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Euler Parameter $\epsilon_4$","Interpreter","latex")
title("$\epsilon_4$(t)","Interpreter","latex")

%% Part d1

% error metric
delta = vecnorm(epsilon,2,2) - 1; % [-]

% plot
figure(2)
plot(t,delta,'r')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Error Metric $\delta$","Interpreter","latex")
title("$\delta(t)$","Interpreter","latex")

%% Part d2 (1e-4)

% new integration options
opts = odeset("RelTol",1e-4,"AbsTol",1e-10);

% integration
[t,epsilon] = ode45(@odefun1,tspan,epsilon0,opts);

% new error metric
delta = vecnorm(epsilon,2,2) - 1; % [-]

% new plot
figure(3)
plot(t,delta,'r')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Error Metric $\delta$","Interpreter","latex")
title("$\delta(t)$ with tolerance of $1\times10^{-4}$","Interpreter","latex")

%% Part d2 (1e-6)

% new integration options
opts = odeset("RelTol",1e-6,"AbsTol",1e-10);

% integration
[t,epsilon] = ode45(@odefun1,tspan,epsilon0,opts);

% new error metric
delta = vecnorm(epsilon,2,2) - 1; % [-]

% new plot
figure(4)
plot(t,delta,'r')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Error Metric $\delta$","Interpreter","latex")
title("$\delta(t)$ with tolerance of $1\times10^{-6}$","Interpreter","latex")

%% Part d3 (1e-8)

% new integration options
opts = odeset("RelTol",1e-8,"AbsTol",1e-10);

% integration
[t,epsilon] = ode45(@odefun1,tspan,epsilon0,opts);

% new error metric
delta = vecnorm(epsilon,2,2) - 1; % [-]

% new plot
figure(5)
plot(t,delta,'r')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Error Metric $\delta$","Interpreter","latex")
title("$\delta(t)$ with tolerance of $1\times10^{-8}$","Interpreter","latex")

%% Problem 02

clearvars -except t epsilon
clc

%% Part b

% initialize variables
theta1 = zeros(size(t)); % [deg]
theta2 = zeros(size(t)); % [deg]
theta3 = zeros(size(t)); % [deg]

% loop over all the EP
for i = 1:length(t)

    % convert EP to DCM
    C = EPtoDCM(epsilon(i,:)); % [-]

    % convert DCM to (3-1-3) sequence Euler Angles
    [alpha,beta,gamma] = DCMtoEA313(C); % [rad]

    % convert angles to deg over interval [-180,180]
    theta1(i) = mod(rad2deg(alpha),360) - 180; % [deg]
    theta2(i) = mod(rad2deg(beta),360) - 180; % [deg]
    theta3(i) = mod(rad2deg(gamma),360) - 180; % [deg]
    
end

% plots
figure(6)
subplot(3,1,1)
plot(t,theta1,'b')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Euler Angle $\theta_1$ [deg]","Interpreter","latex")
title("$\theta_1(t)$","Interpreter","latex")
ylim([-180 180])

subplot(3,1,2)
plot(t,theta2,'b')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Euler Angle $\theta_2$ [deg]","Interpreter","latex")
title("$\theta_2(t)$","Interpreter","latex")
ylim([-180 180])

subplot(3,1,3)
plot(t,theta3,'b')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Euler Angle $\theta_3$ [deg]","Interpreter","latex")
title("$\theta_3(t)$","Interpreter","latex")
ylim([-180 180])


%% Part d

% time span
tspan = 0:0.5:100; % [s]

% integration options
opts = odeset("RelTol",1e-10,"AbsTol",1e-10);

% initial theta vector
theta0 = [-pi/4; pi/8; pi/5]; % [rad]

% integration
[t,theta] = ode45(@odefun2,tspan,theta0,opts);

% convert to degrees over [-180, 180]
theta = mod(rad2deg(theta), 360) - 180; % [deg]

% plots
figure(7)
subplot(3,1,1)
plot(t,theta(:,1),'b')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Euler Angle $\theta_1$ [deg]","Interpreter","latex")
title("$\theta_1(t)$","Interpreter","latex")

subplot(3,1,2)
plot(t,theta(:,2),'b')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Euler Angle $\theta_2$ [deg]","Interpreter","latex")
title("$\theta_2(t)$","Interpreter","latex")

subplot(3,1,3)
plot(t,theta(:,3),'b')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Euler Angle $\theta_3$ [deg]","Interpreter","latex")
title("$\theta_3$(t)","Interpreter","latex")


%% part e

theta_dot = zeros(size(theta));

for i = 1:length(t)
    theta_dot(i,:) = odefun2(t(i),theta(i,:));
end % for

% plots
figure(8)
subplot(3,1,1)
plot(t,theta_dot(:,1),'b')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Euler Angle $\dot{\theta_1}$ [deg/s]","Interpreter","latex")
title("$\dot{\theta_1}(t)$","Interpreter","latex")

subplot(3,1,2)
plot(t,theta_dot(:,2),'b')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Euler Angle $\dot{\theta_2}$ [deg/s]","Interpreter","latex")
title("$\dot{\theta_2}(t)$","Interpreter","latex")

subplot(3,1,3)
plot(t,theta_dot(:,3),'b')
grid on
xlabel("Time [s]","Interpreter","latex")
ylabel("Euler Angle $\dot{\theta_3}$ [deg/s]","Interpreter","latex")
title("$\dot{\theta_3}$(t)","Interpreter","latex")

%% ODE for part c (problem 1)
function epsilon_dot = odefun1(t,epsilon)

% current angular rate
omega(1) = 0.1*cos(0.2*t) + 0.15*sin(0.2*t); % [rad/s]
omega(2) = 0.15*cos(0.2*t) - 0.1*sin(0.2*t); % [rad/s]
omega(3) = 0.3; % [rad/s]


% Euler Parameter rate
epsilon_dot = KDE_EP(epsilon,omega); % [-]

end

%% ODE for part d (problem 2)
function theta_dot = odefun2(t,theta)

% current angular rate
omega = zeros(3,1);
omega(1) = 0.1*cos(0.2*t) + 0.15*sin(0.2*t); % [rad/s]
omega(2) = 0.15*cos(0.2*t) - 0.1*sin(0.2*t); % [rad/s]
omega(3) = 0.3; % [rad/s]

% Euler angle rate
theta_dot = KDE_EA313(theta,omega);

end
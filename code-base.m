clc;
clear;
close all;

% Van der Pol system:
% x' = y
% y' = mu*(1 - x^2)*y - x

mu = 2;

% Initial conditions
t0 = 1;
x0 = 3;
y0 = 5;

% Final time
T = 3;

% Step for Euler method
h_euler = 0.025;

% Step for Runge-Kutta method
h_rk = 0.1;

%% ============================
% Euler Method
%% ============================

t_euler = t0:h_euler:T;
N_euler = length(t_euler);

x_euler = zeros(1, N_euler);
y_euler = zeros(1, N_euler);

x_euler(1) = x0;
y_euler(1) = y0;

for n = 1:N_euler-1
    f = y_euler(n);
    g = mu * (1 - x_euler(n)^2) * y_euler(n) - x_euler(n);

    x_euler(n+1) = x_euler(n) + h_euler * f;
    y_euler(n+1) = y_euler(n) + h_euler * g;
end

%% ============================
% Runge-Kutta 4th Order Method
%% ============================

t_rk = t0:h_rk:T;
N_rk = length(t_rk);

x_rk = zeros(1, N_rk);
y_rk = zeros(1, N_rk);

x_rk(1) = x0;
y_rk(1) = y0;

for n = 1:N_rk-1
    t = t_rk(n);
    x = x_rk(n);
    y = y_rk(n);

    % k1
    k1x = y;
    k1y = mu * (1 - x^2) * y - x;

    % k2
    x2 = x + h_rk/2 * k1x;
    y2 = y + h_rk/2 * k1y;

    k2x = y2;
    k2y = mu * (1 - x2^2) * y2 - x2;

    % k3
    x3 = x + h_rk/2 * k2x;
    y3 = y + h_rk/2 * k2y;

    k3x = y3;
    k3y = mu * (1 - x3^2) * y3 - x3;

    % k4
    x4 = x + h_rk * k3x;
    y4 = y + h_rk * k3y;

    k4x = y4;
    k4y = mu * (1 - x4^2) * y4 - x4;

    % Update
    x_rk(n+1) = x + h_rk/6 * (k1x + 2*k2x + 2*k3x + k4x);
    y_rk(n+1) = y + h_rk/6 * (k1y + 2*k2y + 2*k3y + k4y);
end

%% ============================
% Display Results
%% ============================

fprintf('Euler method result at t = %.2f:\n', T);
fprintf('x(%.2f) = %.4f\n', T, x_euler(end));
fprintf('y(%.2f) = %.4f\n\n', T, y_euler(end));

fprintf('Runge-Kutta 4th order result at t = %.2f:\n', T);
fprintf('x(%.2f) = %.4f\n', T, x_rk(end));
fprintf('y(%.2f) = %.4f\n', T, y_rk(end));

%% ============================
% Plot x(t)
%% ============================

figure;
plot(t_euler, x_euler, 'LineWidth', 1.5);
hold on;
plot(t_rk, x_rk, 'LineWidth', 1.5);
grid on;

xlabel('t');
ylabel('x(t)');
title('Euler and Runge-Kutta Methods for Van der Pol Oscillator');
legend('Euler Method', 'Runge-Kutta 4th Order');

%% ============================
% Phase Portrait
%% ============================

figure;
plot(x_euler, y_euler, 'LineWidth', 1.5);
hold on;
plot(x_rk, y_rk, 'LineWidth', 1.5);
grid on;

xlabel('x');
ylabel('y');
title('Phase Portrait');
legend('Euler Method', 'Runge-Kutta 4th Order');

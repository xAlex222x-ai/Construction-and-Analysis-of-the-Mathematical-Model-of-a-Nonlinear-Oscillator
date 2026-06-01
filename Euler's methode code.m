clc;
clear;
close all;

mu = 2;

t0 = 1;
T = 3;

x0 = 3;
y0 = 5;

eps_required = 0.1;

h = 0.4;

max_iter = 20;

fprintf('Euler method with automatic step halving:\n');
fprintf(' h        x_h(3)          y_h(3)          x_h/2(3)        y_h/2(3)        R\n');

for iter = 1:max_iter

    [x_h, y_h, t_values_h, x_values_h, y_values_h] = EulerSolve(mu,t0,T,x0,y0,h);

    h_half = h/2;

    [x_half, y_half, t_values_half, x_values_half, y_values_half] = EulerSolve(mu,t0,T,x0,y0,h_half);

    Rx = abs(x_half - x_h);
    Ry = abs(y_half - y_h);

    R = max(Rx,Ry);

    fprintf('%6.4f   %12.4e   %12.4e   %12.4e   %12.4e   %12.4e\n', ...
        h, x_h, y_h, x_half, y_half, R);

    if R < eps_required

        fprintf('\nAccuracy condition is satisfied: R = %.4f < eps = %.4f\n', R, eps_required);
        fprintf('The step is accepted: h = %.4f\n', h_half);
        fprintf('Final result:\n');
        fprintf('x(3) = %.4f\n', x_half);
        fprintf('y(3) = %.4f\n', y_half);

        t_euler = t_values_half;
        x_euler = x_values_half;
        y_euler = y_values_half;

        break;

    else

        h = h_half;

    end

end

figure;
plot(t_euler,x_euler,'LineWidth',1.5);
grid on;
xlabel('t');
ylabel('x(t)');
title('Euler Method for Van der Pol Oscillator');

figure;
plot(x_euler,y_euler,'LineWidth',1.5);
grid on;
xlabel('x');
ylabel('y');
title('Phase Portrait, Euler Method');

function [x_end, y_end, t, x, y] = EulerSolve(mu,t0,T,x0,y0,h)

    t = t0:h:T;
    N = length(t);

    x = zeros(1,N);
    y = zeros(1,N);

    x(1) = x0;
    y(1) = y0;

    for n = 1:N-1

        f = y(n);
        g = mu*(1 - x(n)^2)*y(n) - x(n);

        x(n+1) = x(n) + h*f;
        y(n+1) = y(n) + h*g;

    end

    x_end = x(end);
    y_end = y(end);

end

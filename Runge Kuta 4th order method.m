clc;
clear;
close all;

mu = 2;

t0 = 1;
T = 1.4;

x0 = 3;
y0 = 5;

h = 0.1;

t = t0:h:T;
N = length(t);

x = zeros(1,N);
y = zeros(1,N);

x(1) = x0;
y(1) = y0;

for n = 1:N-1

    Xn = [x(n); y(n)];

    k1 = F(Xn,mu);
    k2 = F(Xn + h/2*k1,mu);
    k3 = F(Xn + h/2*k2,mu);
    k4 = F(Xn + h*k3,mu);

    X_next = Xn + h/6*(k1 + 2*k2 + 2*k3 + k4);

    x(n+1) = X_next(1);
    y(n+1) = X_next(2);

end

fprintf('Runge-Kutta 4th order method:\n');
fprintf(' n      t        x(t)        y(t)\n');

for i = 1:N
    fprintf('%2d   %4.1f   %10.4f   %10.4f\n', ...
        i-1, t(i), x(i), y(i));
end

fprintf('\nFinal result:\n');
fprintf('x(%.1f) = %.4f\n',T,x(end));
fprintf('y(%.1f) = %.4f\n',T,y(end));

figure;
plot(t,x,'o-','LineWidth',1.5);
grid on;
xlabel('t');
ylabel('x(t)');
title('Runge-Kutta 4th Order Method for Van der Pol Oscillator');

figure;
plot(x,y,'o-','LineWidth',1.5);
grid on;
xlabel('x');
ylabel('y');
title('Phase Portrait, Runge-Kutta 4th Order Method');

function dX = F(X,mu)

    x = X(1);
    y = X(2);

    dX = [
        y;
        mu*(1 - x^2)*y - x
    ];

end

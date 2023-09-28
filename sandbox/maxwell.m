
clear
clc

X = 0.1;
freq = 10.0;
omega = 2 * pi * freq;
nrm = 100.0;

t = linspace(0.0, 1.0, 1000);
x = X * cos(omega * t);
v = -X * omega * sin(omega * t);

mdl.friction_coefficient = 0.3;
mdl.stiffness = 1.1e3;

F = zeros(size(t));
xi = 0.0;
di = 0.0;
delta = mdl.friction_coefficient * nrm / mdl.stiffness
for i = 1:length(t)
    zeta = x(i) - xi + di;
    di = sign(zeta) * min(abs(zeta), delta);
    F(i) = mdl.stiffness * di;
    xi = x(i);
end

subplot(3, 1, 1)
plot(t, F)
ylabel('F(t)')

subplot(3, 1, 2)
plot(t, x)
ylabel('x(t)')

subplot(3, 1, 3)
plot(t, v)
ylabel('v(t)')
pause

subplot(1, 1, 1)
plot(x, F)
xlabel('x(t)')
ylabel('F(t)')
pause

plot(v, F)
xlabel('v(t)')
ylabel('F(t)')
pause
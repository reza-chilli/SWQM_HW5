clc; clear;

depth = 4; % m
width_es = 100; % m
velocity = 250; % m/hr
E = 1000000; % m^2/hr
W = 10 * 1000; % gr
k = 0.1/24; % 1/hr
area = depth * width_es; % m^2

W_plane = W/area; % gr/m^2

deltaT = 1; % hr
deltaX = 100; % m

t = 1:deltaT:10; % hr (inspect in a day with interval of 3 hr)
x = 0:deltaX:10000; % m (inspect in 10 km every 20 m)

% applying analytical method
c_analytical = zeros(length(t), length(x)); % gr/m^3 or ppm

for i = 1:size(c_analytical, 1)
    for j = 1:size(c_analytical, 2)
        % eq 10.24 chapra:
        c_analytical(i, j) = (W_plane/(2 * ((pi * E * t(i))^ 0.5))) * exp(-(((x(j) - ((velocity * t(i)))) ^ 2)/(4 * E * t(i))) - (k * t(i)));
    end
end

% initiating plot with labels
figure;
plot(x, c_analytical(1, :));
hold on;
plot(x, c_analytical(2, :));
plot(x, c_analytical(4, :));
plot(x, c_analytical(6, :));
plot(x, c_analytical(10, :));

xlabel('Length (m)');
ylabel('Concentration (ppm)');
legend('t = 1h', 't = 2h', 't = 4h', 't = 6h', 't = 10h');
hold off;
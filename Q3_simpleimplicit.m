clc; clear;

depth = 4; % m
width_es = 100; % m
velocity = 250; % m/hr
E = 1000000; % m^2/hr
W = 10 * 1000; % gr
k = 0.1/24; % 1/hr
area = depth * width_es; % m^2

% simple implicit method (applying centered difference with alpha = 0.5 and
% beta = 0.5
deltaT = 0.1; % 0.1 hour
deltaX = 500; % 500m
x = 0:deltaX:10000;
t = 0:deltaT:24; % 0.1h time step

c_implicit = zeros(length(t) + 1, length(x)); % gr/m^3 or ppm
c_implicit(1, 1) = W/(area * deltaX); % setting initial conditions

E = E - (((velocity ^ 2) * deltaT)/2); % applying numerical dispersion
coefficients = zeros(size(c_implicit, 2), size(c_implicit, 2));

for j = 1:size(c_implicit, 2)
    coefficients(j, j) = 1 + ((2 * E * deltaT)/(deltaX ^ 2)) + (k * deltaT);
    switch j
        case 1
            coefficients(j, j + 1) = ((velocity * deltaT)/(2 * deltaX)) - ((E * deltaT)/(deltaX ^ 2));
        case size(c_implicit, 2)
            coefficients(j, j - 1) = -((velocity * deltaT)/(2 * deltaX)) - ((E * deltaT)/(deltaX ^ 2));
        otherwise
            coefficients(j, j - 1) = -((velocity * deltaT)/(2 * deltaX)) - ((E * deltaT)/(deltaX ^ 2));
            coefficients(j, j + 1) = ((velocity * deltaT)/(2 * deltaX)) - ((E * deltaT)/(deltaX ^ 2));
    end
end

for i = 2:size(c_implicit, 1)
    c_implicit(i, :) = coefficients\flipud(rot90(c_implicit(i - 1, :)));
end

% initiating plot with labels
figure;
plot(x, c_implicit(11, :));
hold on;
plot(x, c_implicit(21, :));
plot(x, c_implicit(41, :));
plot(x, c_implicit(61, :));
plot(x, c_implicit(101, :));

xlabel('Length (m)');
ylabel('Concentration (ppm)');
legend('t = 1h', 't = 2h', 't = 4h', 't = 6h', 't = 10h');
hold off;
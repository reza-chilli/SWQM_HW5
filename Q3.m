clc; clear;

depth = 4; % m
width_es = 100; % m
velocity = 250/(60 * 60); % m/s
E = 1000000/(60 * 60); % m^2/s
W = 10; % kg
k = 0.1/(24 * 60 * 60); % 1/s
area = depth * width_es; % m^2
discharge = width_es * depth * velocity/(60 * 60); % m^3/s

W_plane = (W/area) * 1000; % gr/m^2

t = 10800:25200:86400; % sec (inspect in a day with interval of 3 hr)
x = 0:20:10000; % m (inspect in 10 km every 20 m)

c_analytical = zeros(length(t), length(x)); % gr/m^3 or ppm

for i = 1:size(c_analytical, 1)
    for j = 1:size(c_analytical, 2)
        c_analytical(i, j) = (W_plane/(2 * ((pi * E * t(i))^ 0.5))) * exp(-(((x(j) - ((velocity * t(i)))) ^ 2)/(4 * E * t(i))) - (k * t(i)));
    end
end

% initiating plot with labels
figure;
plot(x, c_analytical(1, :));
hold on;

for i = 2:size(c_analytical, 1)
    plot(x, c_analytical(i, :));
end
xlabel('Length (m)');
ylabel('Concentration (ppm)');
legend('t = 3h', 't = 10h', 't = 17h', 't = 24h');
hold off;



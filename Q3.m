clc; clear;

depth = 4; % m
width_es = 100; % m
velocity = 250; % m/hr
E = 1000000; % m^2/hr
W = 10; % kg
k = 0.1/24; % 1/hr
area = depth * width_es; % m^2
discharge = width_es * depth * velocity; % m^3/hr

W_plane = W/area; % kg/m^2

t = 3:7:24; % hr (inspect in a day with interval of 3 hr)
x = 0:0.25:10; % m (inspect in 1 km every 20 m)

c_analytical = zeros(length(t), length(x));

for i = 1:size(c_analytical, 1)
    for j = 1:size(c_analytical, 2)
        c_analytical(i, j) = (W_plane/(2 * ((pi * 0.1 * t(i))^ 0.5))) * exp(-(((x(j) - ((0.2 * t(i)))) ^ 2)/(4 * 0.1 * t(i))));
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
ylabel('Concentration (kg/m^3)');
legend('t = 3', 't = 10', 't = 17', 't = 24');
hold off;



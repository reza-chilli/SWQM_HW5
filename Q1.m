width = 2; % m
height = 1; % m
area = width * height; % m^2
length = 10; % m
k = 0.2; % 1/yr
E = 2; % m^2/hr
c_in = 100; % mg/L
velocity = 1; % m/hr
% eq 9.25 chapra
eta = (k * E)/(velocity ^ 2);
lambda_1 = (velocity/(2 * E)) * (1 + ((1 + (4 * eta))^0.5));
lambda_2 = (velocity/(2 * E)) * (1 - ((1 + (4 * eta))^0.5));
% eq 9.31 chapra
F = (velocity * c_in * lambda_2 * (exp(lambda_2 * length)))/((velocity - (E * lambda_1)) *  lambda_2 * (exp(lambda_2 * length)) - ((velocity - (E * lambda_2)) * lambda_1 * (exp(lambda_1 * length))));
% eq 9.32 chapra
G = (velocity * c_in * lambda_1 * (exp(lambda_1 * length)))/((velocity - (E * lambda_2)) *  lambda_1 * (exp(lambda_1 * length)) - ((velocity - (E * lambda_1)) * lambda_2 * (exp(lambda_2 * length))));
x = 0:0.5:length;
% eq 9.26 chapra
analytical = (F * exp(lambda_1 * x)) + (G * exp(lambda_2 * x));
plot(x, analytical);
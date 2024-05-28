width = 2; % m
height = 1; % m
area = width * height; % m^2
L = 10; % m
k = 0.2; % 1/hr
E = 2; % m^2/hr
c_in = 100; % mg/L
velocity = 1; % m/hr
% eq 9.25 chapra
eta = (k * E)/(velocity ^ 2);
lambda_1 = (velocity/(2 * E)) * (1 + ((1 + (4 * eta))^0.5));
lambda_2 = (velocity/(2 * E)) * (1 - ((1 + (4 * eta))^0.5));
% eq 9.31 chapra
F = (velocity * c_in * lambda_2 * (exp(lambda_2 * L)))/((velocity - (E * lambda_1)) *  lambda_2 * (exp(lambda_2 * L)) - ((velocity - (E * lambda_2)) * lambda_1 * (exp(lambda_1 * L))));
% eq 9.32 chapra
G = (velocity * c_in * lambda_1 * (exp(lambda_1 * L)))/((velocity - (E * lambda_2)) *  lambda_1 * (exp(lambda_1 * L)) - ((velocity - (E * lambda_1)) * lambda_2 * (exp(lambda_2 * L))));
x = 0:0.5:L; % assuming n = 20
% eq 9.26 chapra
analyticalSolution = (F * exp(lambda_1 * x)) + (G * exp(lambda_2 * x));
% backward solution
n = length(x) - 1;
coefficients = zeros(n, n);
loads = zeros(n, 1);

E_prime = (E * area)/(L/n); % E_prime is the same for all elements no need to calculate in loop
discharge = velocity * area; % discharge is the same for all elements no need to calculate in loop
elementVolume = area * (L/n); % element volume is the same for all elements no need to calculate in loop
loads(1, 1) = c_in * discharge;
for i = 1:1:n
    switch i
        case 1
            coefficients(i, i) = discharge + E_prime + (k * elementVolume);
            coefficients(i, i + 1) = -E_prime;
        case n
            coefficients(i, i - 1) = -discharge - E_prime;
            coefficients(i, i) = discharge + E_prime + (k * elementVolume);
        otherwise
            coefficients(i, i - 1) = -discharge - E_prime;
            coefficients(i, i) = discharge + (2 * E_prime) + (k * elementVolume);
            coefficients(i, i + 1) = -E_prime;
    end
end

backwardSolution = coefficients\loads;
% initiating plot with labels
figure;
plot(x, analyticalSolution);
hold on;

centered_x = 0.25:0.5:9.75;
plot(centered_x, backwardSolution, '-o');

xlabel('Length (m)');
ylabel('Concentration (mg/L)');
legend('Analytical Solution', 'Backward Solution');
hold off;
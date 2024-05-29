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
% backward and centered solution
n = length(x) - 1;
backwardCoefficients = zeros(n, n);
centeredCoefficients = zeros(n, n);
loads = zeros(n, 1);

E_prime = (E * area)/(L/n); % E_prime is the same for all elements no need to calculate in loop
discharge = velocity * area; % discharge is the same for all elements no need to calculate in loop
elementVolume = area * (L/n); % element volume is the same for all elements no need to calculate in loop
loads(1, 1) = c_in * discharge;
for i = 1:1:n
    switch i
        case 1
            backwardCoefficients(i, i) = discharge + E_prime + (k * elementVolume);
            centeredCoefficients(i, i) = E_prime + (k * elementVolume) + (discharge/2);
            backwardCoefficients(i, i + 1) = -E_prime;
            centeredCoefficients(i, i + 1) = -E_prime + (discharge/2);
        case n
            backwardCoefficients(i, i - 1) = -discharge - E_prime;
            centeredCoefficients(i, i - 1) = -(discharge/2) - E_prime;
            backwardCoefficients(i, i) = discharge + E_prime + (k * elementVolume);
            centeredCoefficients(i, i) = (2 * E_prime) + (k * elementVolume) - (discharge/2);
        otherwise
            backwardCoefficients(i, i - 1) = -discharge - E_prime;
            centeredCoefficients(i, i - 1) = -(discharge/2) - E_prime;
            backwardCoefficients(i, i) = discharge + (2 * E_prime) + (k * elementVolume);
            centeredCoefficients(i, i) = (2 * E_prime) + (k * elementVolume);
            backwardCoefficients(i, i + 1) = -E_prime;
            centeredCoefficients(i, i + 1) = -E_prime + (discharge/2);
    end
end

backwardSolution = backwardCoefficients\loads;
centeredSolution = centeredCoefficients\loads;
% initiating plot with labels
figure;
plot(x, analyticalSolution);
hold on;

centered_x = 0.25:0.5:9.75;
plot(centered_x, backwardSolution, '-o');
plot(centered_x, centeredSolution, '-x');
xlabel('Length (m)');
ylabel('Concentration (mg/L)');
legend('Analytical Solution', 'Backward Solution', 'Centered Solution');
hold off;
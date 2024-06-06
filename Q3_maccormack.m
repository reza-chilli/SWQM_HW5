clc; clear;

depth = 4; % m
width_es = 100; % m
velocity = 250; % m/hr
E = 1000000; % m^2/hr
W = 10 * 1000; % gr
k = 0.1/24; % 1/hr
area = depth * width_es; % m^2

% applying maccormack method
deltaT = 0.00001; % hr
deltaX = 100; % m
t = 0:deltaT:10;
x = 0:deltaX:10000;
% setting concentration matrix, each row is associated to a time step and
% each column is associated to an element
% column 1, n are associated to boundry conditions
% row 1 is associated to initial conditions
c_maccormack = zeros(length(t) + 1, length(x) + 2); % gr/m^3 or ppm (index 1, n is associated to boundry conditions)
c_maccormack(1, 2) = W/(area * deltaX); % setting initial conditions

for i = 2:size(c_maccormack, 1)
    % slopes matrix first row is associated to predictor, second row is
    % associated to corrector
    slopes = zeros(2, size(c_maccormack, 2));
    dummyConcentrations = zeros(1, size(c_maccormack, 2) + 2);
    for j = 2:size(c_maccormack, 2) - 1
        % eq 13.12 chapra: 
        slopes(1, j) = -velocity * (c_maccormack(i-1, j+1) - c_maccormack(i-1, j)) / deltaX + E * (c_maccormack(i-1, j+1) - 2 * c_maccormack(i-1, j) + c_maccormack(i-1, j-1)) / (deltaX^2) - k * c_maccormack(i-1, j);
        % eq 13.13 chapra:
        dummyConcentrations(1, j) = c_maccormack(i - 1, j) + (slopes(1, j) * deltaT);
    end
    for j = 2:size(c_maccormack, 2) - 1
        % eq 13.14 chapra
        slopes(2, j) = -velocity * (dummyConcentrations(j) - dummyConcentrations(j-1)) / deltaX + E * (dummyConcentrations(j+1) - 2 * dummyConcentrations(j) + dummyConcentrations(j-1)) / (deltaX^2) - k * dummyConcentrations(j);
        % eq 13.15 chapra
        c_maccormack(i, j) = c_maccormack(i - 1, j) + (((slopes(1, j) + slopes(2, j))/2) * deltaT);
    end
end
figure;
plot(x, c_maccormack(100001, 2:size(c_maccormack, 2) - 1)); % 1hr
hold on;
plot(x, c_maccormack(200001, 2:size(c_maccormack, 2) - 1)); % 2hr
plot(x, c_maccormack(400001, 2:size(c_maccormack, 2) - 1)); % 4hr
plot(x, c_maccormack(600001, 2:size(c_maccormack, 2) - 1)); % 6hr
plot(x, c_maccormack(1000001, 2:size(c_maccormack, 2) - 1)); % 10hr
xlabel('Length (m)');
ylabel('Concentration (ppm)');
legend('t = 1h', 't = 2h', 't = 4h', 't = 6h', 't = 10h');
hold off;

clc; clear;

A_surface = [280, 280, 280]; % km^2
depth = [4.3, 4.3, 4.3]; % m

segmentVolume = zeros(1, 3); % km^3
for i = 1:length(segmentVolume)
    segmentVolume(i) = A_surface(i) * depth(i)/1000;
end

W = 20; % kg into segment 2 (since one load is applied no need for an array)
W_s = 0.02; % km/yr
Q_12 = 170; % km^3/yr
Q_23 = 10; % km^3/yr
Q_2_out = 180; % km^3/yr
E_prime_12 = 6.5; % km^3/yr
E_prime_23 = 6.5; % km^3/yr

length_12 = 16.7; % km
length_23 = 16.7; % km
A_cross_12 = length_12 * ((depth(1) + depth(2))/2000); % km^2
A_cross_23 = length_23 * ((depth(2) + depth(3))/2000); % km^2

timeStep = 0.05; % day
inspectedTime = 20; % day

c_1 = zeros(1, inspectedTime/timeStep); % concentration of segment 1 over time (kg/km^3)
c_2 = zeros(1, inspectedTime/timeStep); % concentration of segment 2 over time (kg/km^3)
c_3 = zeros(1, inspectedTime/timeStep); % concentration of segment 3 over time (kg/km^3)

c_2(1) = W/(segmentVolume(2)); % kg/km^3

% apply explicit FTBS method
alpha = 1;
beta = 0;

for i = 2:(inspectedTime/timeStep) + 1
    % eq 12.13 chapra (modified)
    c_1(i) = ((1 - ((timeStep/(segmentVolume(1) * 365)) * ((Q_12 * alpha) + E_prime_12 + (W_s * A_surface(1))))) * c_1(i - 1)) - ((timeStep/(segmentVolume(1) * 365)) * (((Q_12 * beta) - E_prime_12) * c_2(i - 1)));
    c_2(i) = -((timeStep/(segmentVolume(2) * 365)) * (-(Q_12 * alpha) - E_prime_12) * c_1(i - 1)) + ((1 - ((timeStep/(segmentVolume(2) * 365)) * (-(Q_12 * beta) + (Q_23 * alpha) + E_prime_12 + E_prime_23 + (W_s * A_surface(2)) + Q_2_out))) * c_2(i - 1)) - ((timeStep/(segmentVolume(2) * 365)) * (((Q_23 * beta) - E_prime_23) * c_3(i - 1)));
    c_3(i) = -((timeStep/(segmentVolume(3) * 365)) * (-(Q_23 * alpha) - E_prime_23) * c_2(i - 1)) + ((1 - ((timeStep/(segmentVolume(3) * 365)) * (-(Q_23 * beta) + E_prime_23 + (W_s * A_surface(3))))) * c_3(i - 1));
end

% initiating plot with labels
figure;
x = 0:timeStep:inspectedTime;
plot(x, c_1);
hold on;

plot(x, c_2);
plot(x, c_3);
xlabel('Time (Day)');
ylabel('Concentration (kg/km^3)');
legend('Segment 1', 'Segment 2', 'Segment 3');
hold off;



% Clear the workspace, close figures, and clear the command window
clear; close all; clc;

% Load the dataset from the specified .mat file
data = load('Candidate_249383.mat');

% Extract variables from the loaded data
Time_Minutes = data.Time_Minutes;
Supply_Flows = data.Supply_Flows;

% Verify the dimensions of the data
assert(numel(Time_Minutes) == 1440, 'Time_Minutes must contain 1440 elements.');
assert(size(Supply_Flows, 2) == 1440, 'Supply_Flows must have 1440 columns.');

% Define the range for time in minutes
timeRange = 1:1440;

% Gaussian fit coefficients for different datasets
fitCoefficients = [
    23.7376, 672.7400, 69.2472;
    56.3213, 815.9610, 134.0388;
    178.9957, 530.8720, 252.9244;
    38.5162, 792.5321, 103.2703;
    149.3198, 544.0592, 239.5997;
    39.0584, 792.5692, 101.4231
];

% Create Gaussian function handles for each dataset
gaussianFuncs = arrayfun(@(i) ...
    @(t) fitCoefficients(i, 1) * exp(-((t - fitCoefficients(i, 2)) / fitCoefficients(i, 3)).^2), ...
    1:6, 'UniformOutput', false);

% Plot the noisy data and Gaussian fits
figure; hold on;
for i = 1:6
    plot(Time_Minutes, Supply_Flows(i, :), '--', 'DisplayName', sprintf('Supply %d (m^3/min)', i));
    plot(timeRange, gaussianFuncs{i}(timeRange), '-', 'DisplayName', sprintf('Supply %d fit', i));
end

% Label the axes and title the plot
xlabel('Time (Minutes)');
ylabel('Supply Flows (m^3/min)');
title('Curve Fitting of Noisy Data');
legend('show');
hold off;

% Define offsets and create adjusted Gaussian function handles
offsets = [-20, -50, -10, -20, -90, -80];
adjustedFuncs = arrayfun(@(i) ...
    @(t) fitCoefficients(i, 1) * exp(-((t - fitCoefficients(i, 2)) / fitCoefficients(i, 3)).^2) - offsets(i), ...
    1:6, 'UniformOutput', false);

% Find the roots of the adjusted Gaussian functions
rootVals = zeros(6, 2);
for i = 1:6
    try
        [a, b] = findInterval(adjustedFuncs{i}, 1, fitCoefficients(i, 2), 1000);
        rootVals(i, 1) = ridder(adjustedFuncs{i}, a, b, 1e-6);
    catch
        rootVals(i, 1) = NaN;
    end
    
    try
        [a, b] = findInterval(adjustedFuncs{i}, fitCoefficients(i, 2), 1440, 1000);
        rootVals(i, 2) = ridder(adjustedFuncs{i}, a, b, 1e-6);
    catch
        rootVals(i, 2) = NaN;
    end
end

% Plot the adjusted Gaussian curves and their roots
figure; hold on;
for i = 1:6
    plot(timeRange, gaussianFuncs{i}(timeRange), 'DisplayName', sprintf('Supply %d', i));
    plot(rootVals(i, 1), adjustedFuncs{i}(rootVals(i, 1)), 'rx', 'HandleVisibility', 'off');
    plot(rootVals(i, 2), adjustedFuncs{i}(rootVals(i, 2)), 'rx', 'HandleVisibility', 'off');
end

% Label the axes and title the plot
xlabel('Time (Minutes)');
ylabel('Supply Flows (m^3/min)');
title('Adjusted Gaussian Curves and Roots');
legend('show');
hold off;

% Simulate head loss over the time range
headLossResults = zeros(1440, 6);
for t = 1:1440
    currentSupply = arrayfun(@(i) gaussianFuncs{i}(t), 1:6);
    [Q, headLoss] = SimulateHeadloss(currentSupply(1), currentSupply(2), currentSupply(3), currentSupply(4), currentSupply(5), currentSupply(6));
    headLossResults(t, :) = headLoss(1:6); % Only consider the first 6 head loss values
end

% Plot the simulated head loss over time
figure; hold on;
plot(timeRange, headLossResults);
xlabel('Time (Minutes)');
ylabel('Head Loss');
title('Head Loss Simulation Over Time');
legend(arrayfun(@(i) sprintf('Headloss %d', i), 1:6, 'UniformOutput', false));
hold off;

% Helper function to find a bracketing interval for a root
function [x1, x2] = findInterval(func, lower, upper, points)
    x = linspace(lower, upper, points);
    f = arrayfun(func, x);
    for k = 1:length(f) - 1
        if sign(f(k)) ~= sign(f(k + 1))
            x1 = x(k);
            x2 = x(k + 1);
            return;
        end
    end
    error('No root found in the interval [%f, %f].', lower, upper);
end

% Ridder's method for root finding
function root = ridder(func, x1, x2, tolerance)
    f1 = func(x1);
    f2 = func(x2);
    while abs(x2 - x1) > tolerance
        x3 = (x1 + x2) / 2;
        f3 = func(x3);
        s = sqrt(f3^2 - f1 * f2);
        if s == 0
            break;
        end
        x4 = x3 + (x3 - x1) * sign(f1 - f2) * f3 / s;
        f4 = func(x4);
        if f4 == 0
            root = x4;
            return;
        end
        if sign(f3) == sign(f4)
            if sign(f1) ~= sign(f4)
                x2 = x4;
                f2 = f4;
            else
                x1 = x4;
                f1 = f4;
            end
        else
            x1 = x3;
            x2 = x4;
            f1 = f3;
            f2 = f4;
        end
    end
    root = (x1 + x2) / 2;
end

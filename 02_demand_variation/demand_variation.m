function demand_variation()
% DEMAND_VARIATION  How supply pressure varies as six units' demand changes over a day.
%   1. Fit a smooth Gaussian curve to each unit's noisy per-minute flow data.
%   2. Use Ridder's method to find when each unit's flow exceeds its declared demand.
%   3. Re-solve the network every minute to get head loss across the day.
close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'lib'));

data = load(fullfile(here, 'demand_data.mat'));
Time_Minutes = data.Time_Minutes;
Supply_Flows = data.Supply_Flows;      % 6 units x 1440 minutes
timeRange = 1:1440;

% Gaussian fit a*exp(-((t-b)/c)^2) for each unit: [a b c]
% (fitted to the noisy data with MATLAB's Curve Fitting Toolbox, gauss1 model)
fitCoefficients = [
     23.7376, 672.7400,  69.2472;
     56.3213, 815.9610, 134.0388;
    178.9957, 530.8720, 252.9244;
     38.5162, 792.5321, 103.2703;
    149.3198, 544.0592, 239.5997;
     39.0584, 792.5692, 101.4231];

gaussianFuncs = arrayfun(@(i) ...
    @(t) fitCoefficients(i, 1) * exp(-((t - fitCoefficients(i, 2)) / fitCoefficients(i, 3)).^2), ...
    1:6, 'UniformOutput', false);

figure; hold on;
for i = 1:6
    plot(Time_Minutes, Supply_Flows(i, :), '--', 'DisplayName', sprintf('Unit %d data', i));
    plot(timeRange, gaussianFuncs{i}(timeRange), '-', 'DisplayName', sprintf('Unit %d fit', i));
end
xlim([0 1440]); xlabel('Time (minutes)'); ylabel('Flow (m^3/min)');
title('Gaussian fits to noisy demand data');
legend('show'); hold off;

% Declared demand for each unit (m^3/min). A unit is over its declared
% demand where fit(t) - demand = 0, so each root marks a crossing.
declaredDemand = [20, 50, 10, 20, 90, 80];
adjustedFuncs = arrayfun(@(i) @(t) gaussianFuncs{i}(t) - declaredDemand(i), ...
    1:6, 'UniformOutput', false);

% Two crossings per unit: one on the rising side of the peak, one on the falling side.
% A unit whose peak never reaches its declared demand returns NaN.
rootVals = NaN(6, 2);
for i = 1:6
    try
        [a, b] = findInterval(adjustedFuncs{i}, 1, fitCoefficients(i, 2), 1000);
        rootVals(i, 1) = ridder(adjustedFuncs{i}, a, b, 1e-6);
    catch
    end
    try
        [a, b] = findInterval(adjustedFuncs{i}, fitCoefficients(i, 2), 1440, 1000);
        rootVals(i, 2) = ridder(adjustedFuncs{i}, a, b, 1e-6);
    catch
    end
end

fprintf('%-6s %16s %16s\n', 'Unit', 'Exceeds from', 'Back below at');
for i = 1:6
    fprintf('%-6d %16s %16s\n', i, minutesToClock(rootVals(i, 1)), minutesToClock(rootVals(i, 2)));
end

figure; hold on;
for i = 1:6
    plot(timeRange, gaussianFuncs{i}(timeRange), 'DisplayName', sprintf('Unit %d', i));
    plot(rootVals(i, :), [1 1] * declaredDemand(i), 'kx', 'MarkerSize', 8, 'HandleVisibility', 'off');
end
xlim([0 1440]); xlabel('Time (minutes)'); ylabel('Flow (m^3/min)');
title('Fitted demand and declared-demand crossings (x)');
legend('show'); hold off;

% Re-solve the network for every minute of the day
headLossResults = zeros(1440, 6);
for t = 1:1440
    currentSupply = arrayfun(@(i) gaussianFuncs{i}(t), 1:6);
    [~, headLoss] = simulate_head_loss(currentSupply(1), currentSupply(2), currentSupply(3), ...
                                       currentSupply(4), currentSupply(5), currentSupply(6));
    headLossResults(t, :) = headLoss(1:6);
end

figure;
plot(timeRange, headLossResults);
xlim([0 1440]); xlabel('Time (minutes)'); ylabel('Head loss (m)');
title('Head loss across the day');
legend(arrayfun(@(i) sprintf('Pipe %d', i), 1:6, 'UniformOutput', false));
end


function [x1, x2] = findInterval(func, lower, upper, points)
% Scan [lower, upper] for the first sign change, to bracket a root for Ridder's method.
    x = linspace(lower, upper, points);
    f = arrayfun(func, x);
    for k = 1:length(f) - 1
        if sign(f(k)) ~= sign(f(k + 1))
            x1 = x(k);
            x2 = x(k + 1);
            return;
        end
    end
    error('No root found in [%f, %f].', lower, upper);
end

function root = ridder(func, x1, x2, tolerance)
% Ridder's method (Ridder, 1979): fits an exponential through the bracket
% to get faster convergence than bisection while keeping the root bracketed.
% Structure follows Kiusalaas, Numerical Methods in Engineering with MATLAB.
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
                x2 = x4; f2 = f4;
            else
                x1 = x4; f1 = f4;
            end
        else
            x1 = x3; x2 = x4; f1 = f3; f2 = f4;
        end
    end
    root = (x1 + x2) / 2;
end

function s = minutesToClock(m)
    if isnan(m)
        s = 'never';
    else
        s = sprintf('%02d:%02d', floor(m / 60), floor(mod(m, 60)));
    end
end

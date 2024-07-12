function Portfolio4_249383
    % Clean up workspace and command window
    close all; % Close all figure windows
    clear; % Clear all variables from workspace
    clc; % Clear command window

    % Initial conditions
    initConc = zeros(1, 12); % Initial concentrations for all chemicals in all tanks
    timeSpan = [0 900]; % Time span for the simulation in seconds
    
    % Solve the system of ODEs
    [time, concentrations] = ode45(@calculateConcentrations, timeSpan, initConc);
    
    % Plot the results
    plotConcentrationData(time, concentrations);
end

% Function defining the rate of change for the ODE system
function dYdt = calculateConcentrations(t, Y)
    % Constants for chemical injections and tank parameters
    alphaInjection = 900; % Injection rate for Alpha
    betaInjection = 60; % Injection rate for Beta
    gammaInjection = 800; % Injection rate for Gamma
    tuningParam = 40; % Tuning parameter for the system
    tankVolume = 300; % Volume of each tank

    % Initialise the rate of change vector
    dYdt = zeros(12, 1);

    % Define the continuity equations for each chemical in each tank

    % Tank 1
    % Rate of change of Alpha in Tank 1
    dYdt(1) = (4*Y(4) - 5*Y(1) + alphaInjection * sin(t + pi/2) + 100) / tankVolume;
    % Rate of change of Beta in Tank 1
    dYdt(2) = (4*Y(5) - 5*Y(2) + betaInjection * sin(t)) / tankVolume;
    % Rate of change of Gamma in Tank 1
    dYdt(3) = (4*Y(6) - 5*Y(3) + gammaInjection * sin(t)) / tankVolume;

    % Tank 2
    % Rate of change of Alpha in Tank 2
    dYdt(4) = (2*Y(4) + Y(7) - 8*Y(4)) / tankVolume;
    % Rate of change of Beta in Tank 2
    dYdt(5) = (2*Y(5) + Y(8) - 8*Y(5) + 700 * sin(t) + tuningParam) / tankVolume;
    % Rate of change of Gamma in Tank 2
    dYdt(6) = (2*Y(6) + Y(9) - 8*Y(6)) / tankVolume;

    % Tank 3
    % Rate of change of Alpha in Tank 3
    dYdt(7) = (2*Y(1) + 2*Y(4) - 5*Y(7)) / tankVolume;
    % Rate of change of Beta in Tank 3
    dYdt(8) = (2*Y(2) + 2*Y(5) - 5*Y(8)) / tankVolume;
    % Rate of change of Gamma in Tank 3
    dYdt(9) = (2*Y(3) + 2*Y(6) - 5*Y(9) + 700 * sin(t + pi/2) + 100) / tankVolume;

    % Tank 4
    % Rate of change of Alpha in Tank 4
    dYdt(10) = (Y(1) + 2*Y(4) + 4*Y(7) - 3*Y(10)) / tankVolume;
    % Rate of change of Beta in Tank 4
    dYdt(11) = (Y(2) + 2*Y(5) + 4*Y(8) - 3*Y(11)) / tankVolume;
    % Rate of change of Gamma in Tank 4
    dYdt(12) = (Y(3) + 2*Y(6) + 4*Y(9) - 3*Y(12)) / tankVolume;
end

% Function to plot the concentrations
function plotConcentrationData(time, concentrations)
    % Create a new figure for the plots
    figure;
    hold on; % Hold on to add multiple plots to the same figure

    % Plot the concentrations for each tank and chemical with different styles
    plot(time, concentrations(:, 1), 'r', 'DisplayName', 'Tank 1: \alpha'); % Tank 1: Alpha
    plot(time, concentrations(:, 2), 'g', 'DisplayName', 'Tank 1: \beta'); % Tank 1: Beta
    plot(time, concentrations(:, 3), 'b', 'DisplayName', 'Tank 1: \gamma'); % Tank 1: Gamma
    plot(time, concentrations(:, 4), 'r--', 'DisplayName', 'Tank 2: \alpha'); % Tank 2: Alpha
    plot(time, concentrations(:, 5), 'g--', 'DisplayName', 'Tank 2: \beta'); % Tank 2: Beta
    plot(time, concentrations(:, 6), 'b--', 'DisplayName', 'Tank 2: \gamma'); % Tank 2: Gamma
    plot(time, concentrations(:, 7), 'r:', 'DisplayName', 'Tank 3: \alpha'); % Tank 3: Alpha
    plot(time, concentrations(:, 8), 'g:', 'DisplayName', 'Tank 3: \beta'); % Tank 3: Beta
    plot(time, concentrations(:, 9), 'b:', 'DisplayName', 'Tank 3: \gamma'); % Tank 3: Gamma
    plot(time, concentrations(:, 10), 'r-.', 'DisplayName', 'Tank 4: \alpha'); % Tank 4: Alpha
    plot(time, concentrations(:, 11), 'g-.', 'DisplayName', 'Tank 4: \beta'); % Tank 4: Beta
    plot(time, concentrations(:, 12), 'b-.', 'DisplayName', 'Tank 4: \gamma'); % Tank 4: Gamma

    % Add labels and title to the plot
    xlabel('Time (seconds)'); % X-axis label
    ylabel('Concentration'); % Y-axis label
    title('Tank Concentrations against Time'); % Plot title

    % Display the legend
    legend('show');

    % Finalise plot
    hold off; % Release the hold on the figure
end

function mixing_tanks
% MIXING_TANKS  Concentrations of three chemicals in four connected tanks.
%   A mass balance on each chemical in each tank gives 12 coupled first-order
%   ODEs (dC/dt = (inflow - outflow + injection) / V), solved with ode45.
%   The time span was extended until every concentration settled into a
%   steady periodic state.
    close all;
    initConc = zeros(1, 12);     % all tanks start clean
    timeSpan = [0 900];          % s

    [time, concentrations] = ode45(@tankBalances, timeSpan, initConc);
    plotConcentrations(time, concentrations);
end

function dYdt = tankBalances(t, Y)
% Y(1:3) = alpha, beta, gamma in tank 1; Y(4:6) tank 2; Y(7:9) tank 3; Y(10:12) tank 4.
    alphaInjection = 900;
    betaInjection  = 60;
    gammaInjection = 800;
    tuningParam    = 40;     % B2, the tunable injection into tank 2
    V = 300;                 % tank volume

    dYdt = zeros(12, 1);

    % Tank 1: fed from tank 2, with oscillating injections of all three chemicals
    dYdt(1) = (4*Y(4) - 5*Y(1) + alphaInjection * sin(t + pi/2) + 100) / V;
    dYdt(2) = (4*Y(5) - 5*Y(2) + betaInjection * sin(t)) / V;
    dYdt(3) = (4*Y(6) - 5*Y(3) + gammaInjection * sin(t)) / V;

    % Tank 2: fed from tank 3, beta injection includes the tuning term
    dYdt(4) = (2*Y(4) + Y(7) - 8*Y(4)) / V;
    dYdt(5) = (2*Y(5) + Y(8) - 8*Y(5) + 700 * sin(t) + tuningParam) / V;
    dYdt(6) = (2*Y(6) + Y(9) - 8*Y(6)) / V;

    % Tank 3: fed from tanks 1 and 2, gamma injection
    dYdt(7) = (2*Y(1) + 2*Y(4) - 5*Y(7)) / V;
    dYdt(8) = (2*Y(2) + 2*Y(5) - 5*Y(8)) / V;
    dYdt(9) = (2*Y(3) + 2*Y(6) - 5*Y(9) + 700 * sin(t + pi/2) + 100) / V;

    % Tank 4: fed from tanks 1, 2 and 3
    dYdt(10) = (Y(1) + 2*Y(4) + 4*Y(7) - 3*Y(10)) / V;
    dYdt(11) = (Y(2) + 2*Y(5) + 4*Y(8) - 3*Y(11)) / V;
    dYdt(12) = (Y(3) + 2*Y(6) + 4*Y(9) - 3*Y(12)) / V;
end

function plotConcentrations(time, C)
% Colour = chemical, line style = tank
    colours = {'r', 'g', 'b'};
    styles  = {'-', '--', ':', '-.'};
    names   = {'\alpha', '\beta', '\gamma'};
    figure; hold on;
    for tank = 1:4
        for chem = 1:3
            col = 3 * (tank - 1) + chem;
            plot(time, C(:, col), [colours{chem} styles{tank}], ...
                'DisplayName', sprintf('Tank %d: %s', tank, names{chem}));
        end
    end
    xlabel('Time (s)'); ylabel('Concentration');
    title('Tank concentrations over time');
    legend('show', 'Location', 'eastoutside');
    hold off;
end

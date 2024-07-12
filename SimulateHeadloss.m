function [Q, head_loss] = SimulateHeadloss(S1, S2, S3, S4, S5, S6)
    % Define the lengths, diameters, and other constants
    Length = [300, 50, 50, 300, 50, 50, 50, 50, 50]; % Lengths of the pipes (in meters)
    Diameter = [0.3, 0.1, 0.1, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1]; % Diameters of the pipes (in meters)
    Density = 1000; % Fluid density (in kg/m^3)
    Gravity = 9.81; % Acceleration due to gravity (in m/s^2)
    U = [1, 1, 1, 1, 1, 1, 1, 1, 1]; % Initial guess for the velocity in each pipe section (in m/s)

    % Initialize variables for the iterative process
    i = 1; % Iteration counter
    tol = 6e-6; % Convergence tolerance
    error = 1; % Initial error

    % Start the iterative process
    while error > tol
        % Calculate Reynolds number, friction factor, and other variables
        Re = (Density .* U .* Diameter) / 0.001002; % Reynolds number
        f = 64 ./ Re; % Friction factor for laminar flow
        k = (8 .* f .* Length) ./ (pi^2 .* Gravity .* Diameter.^6); % Head loss coefficient
        q = U .* (pi .* (Diameter.^2) ./ 4); % Flow rate in each pipe

        % Construct the matrices for the linear system
        B = 2 .* k .* q; % Coefficient matrix B
        Bprime = (1 - 2) .* k .* q .^ 2; % Modified coefficients Bprime

        % Coefficient matrix A representing the node-connection topology
        A = [1 -1 0 0 1 0 0 0 0;
             0 1 1 0 0 -1 0 0 0;
             0 0 -1 1 0 0 1 0 0;
             0 0 0 0 -1 0 0 1 0;
             0 0 0 0 0 1 0 -1 -1;
             0 0 0 0 0 0 -1 0 1;
            -B(1) -B(2) B(3) B(4) 0 0 0 0 0;
             0 B(2) 0 0 B(5) B(6) 0 B(8) 0;
             0 0 -B(3) 0 0 -B(6) -B(7) 0 -B(9)];

        % Source term vector S
        S = [S1;
             S2;
             S3;
             S4;
             S5;
             S6;
             Bprime(1) + Bprime(2) - Bprime(3) - Bprime(4);
            -Bprime(2) - Bprime(5) - Bprime(6) - Bprime(8);
             Bprime(3) + Bprime(6) + Bprime(7) + Bprime(9)];

        % Solve the linear system using Gaussian elimination with pivoting
        Q = gaussPiv(A, S)';

        % Update the error and the flow rates
        error = max(abs(q - Q));
        q = Q;
        i = i + 1;
        U = q ./ (pi .* (Diameter.^2) ./ 4); % Update velocities based on new flow rates
    end

    % Calculate head loss for each pipe
    head_loss = (f .* Length ./ Diameter) .* (U.^2) / (2 * Gravity);
end

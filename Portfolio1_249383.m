function [Q] = Portfolio1_249383()
    clc; % Clear command window
    format short; % Set format to short for numerical output

    % Defined Measurements
    D = [0.30, 0.10, 0.10, 0.30, 0.10, 0.10, 0.10, 0.10, 0.10]; % Pipe Diameters (m)
    L = [300, 50, 50, 300, 50, 50, 50, 50, 50]; % Pipe Lengths (m)
    % Service nodes' flow requirements (m^3/s)
    S = [0.02, 0.04, 0.09, 0.03, 0.08, 0.03];
    % Initial conditions for the while loop
    iter = 1; % Iteration counter
    arrQ = (pi * (D).^2) / 4; % Initial array for flow rates (Q) in each pipe
    diff = ones(1, 9); % Initial difference array for convergence check

    % Iterative while loop for convergence
    while any(diff > 1*10^-9)
        % Calculate coefficients K, B, and B_prime for each pipe
        for I = 1:9
            K(I) = (128 * 1.0016 * 10^-3 * L(I)) / (pi * 9.81 * 1000 * 10^3 * arrQ(I) * (D(I)^5)); % Coefficient K
            B(I) = 2 * K(I) * arrQ(I); % Coefficient B
            B_prime(I) = -K(I) * (arrQ(I)^2); % Coefficient B_prime
        end

        % Define matrix A and vector b for the system of linear equations
        A = [1 -1 0 0 -1 0 0 0 0; % Connectivity matrix
             0 1 1 0 0 -1 0 0 0;
             0 0 -1 1 0 0 -1 0 0;
             0 0 0 0 1 0 0 -1 0;
             0 0 0 0 0 1 0 1 1;
             0 0 0 0 0 0 1 0 -1;
             0 B(2) 0 0 -B(5) B(6) 0 -B(8) 0;
             0 0 -B(3) 0 0 -B(6) B(7) 0 B(9);
            -B(1) -B(2) B(3) B(4) 0 0 0 0 0];

        b = [S(1);
             S(2);
             S(3);
             S(4);
             S(5);
             S(6);
             B_prime(8) + B_prime(5) - B_prime(6) - B_prime(2);
             B_prime(6) + B_prime(3) - B_prime(9) - B_prime(7);
             B_prime(2) + B_prime(1) - B_prime(4) - B_prime(3)];

        iter = iter + 1; % Increment iteration counter
        Q = gaussPiv(A, b); % Solve system of equations using Gaussian elimination
        pre(iter, :) = Q;  % Store the current iteration results
        diff = abs(pre(iter, :) - pre(iter - 1, :));  % Calculate the difference from the previous iteration
    end

    % Calculate head loss for each pipe
    for j = 1:9
        h_loss(j) = abs(K(j)) * (abs(Q(j))^2);
    end

    % Calculate cumulative head loss for each service node
    m_loss = [h_loss(1);
              h_loss(1) + h_loss(2);
              h_loss(4);
              h_loss(1) + h_loss(5);
              h_loss(1) + h_loss(2) + h_loss(6);
              h_loss(4) + h_loss(7)];

    % Display head loss for each pipe
    fprintf("     Head loss\n");
    for i = 1:length(h_loss)
        fprintf("Pipe Q%d: %.4f\n", i, h_loss(i));
    end
    fprintf("\n");

    % Calculate and display pressure loss for each service node
    p_loss = m_loss * 9.804; % Pressure loss in kPa
    fprintf("     Pressure loss\n");
    for i = 1:length(p_loss)
        fprintf("Node S%d: %.4f kPa\n", i, p_loss(i));
    end
    fprintf("\n");

    % Display flow rates for each pipe
    fprintf("     Flow rate\n");
    for i = 1:length(Q)
        fprintf('Pipe Q%d: %.4f L/s\n', i, Q(i));
    end

end




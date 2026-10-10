function Q = water_network()
% WATER_NETWORK  Flow rates and pressure losses in a 9-pipe supply network.
%   Six continuity equations (one per node) and three loop equations are
%   solved together. Head loss is non-linear in Q, so each loop equation is
%   linearised about the previous estimate and the 9x9 system is re-solved
%   until the flows stop changing.
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'lib'));

    D = [0.30, 0.10, 0.10, 0.30, 0.10, 0.10, 0.10, 0.10, 0.10];  % pipe diameters (m)
    L = [300, 50, 50, 300, 50, 50, 50, 50, 50];                   % pipe lengths (m)
    S = [0.02, 0.04, 0.09, 0.03, 0.08, 0.03];                     % unit demands (m^3/s)

    iter = 1;
    arrQ = (pi * D.^2) / 4;     % starting guess for each pipe flow
    diff = ones(1, 9);

    while any(diff > 1e-9)
        for I = 1:9
            K(I) = (128 * 1.0016e-3 * L(I)) / (pi * 9.81 * 1000 * 1e3 * arrQ(I) * D(I)^5);
            B(I) = 2 * K(I) * arrQ(I);          % slope of h = K*Q^2 at the current Q
            B_prime(I) = -K(I) * arrQ(I)^2;     % intercept of that tangent line
        end

        % Rows 1-6: continuity at each node. Rows 7-9: head balance round each loop.
        A = [ 1  -1    0    0   -1    0    0    0    0;
              0   1    1    0    0   -1    0    0    0;
              0   0   -1    1    0    0   -1    0    0;
              0   0    0    0    1    0    0   -1    0;
              0   0    0    0    0    1    0    1    1;
              0   0    0    0    0    0    1    0   -1;
              0  B(2)  0    0  -B(5) B(6)  0  -B(8)  0;
              0   0  -B(3)  0    0  -B(6) B(7)  0   B(9);
            -B(1) -B(2) B(3) B(4)  0    0    0    0    0];

        b = [S(:);
             B_prime(8) + B_prime(5) - B_prime(6) - B_prime(2);
             B_prime(6) + B_prime(3) - B_prime(9) - B_prime(7);
             B_prime(2) + B_prime(1) - B_prime(4) - B_prime(3)];

        iter = iter + 1;
        Q = gaussPiv(A, b);
        pre(iter, :) = Q;
        diff = abs(pre(iter, :) - pre(iter - 1, :));
    end

    h_loss = abs(K) .* abs(Q').^2;

    % Head loss from the source to each unit, following the supply path
    m_loss = [h_loss(1);
              h_loss(1) + h_loss(2);
              h_loss(4);
              h_loss(1) + h_loss(5);
              h_loss(1) + h_loss(2) + h_loss(6);
              h_loss(4) + h_loss(7)];
    p_loss = m_loss * 9.804;    % kPa

    fprintf('Converged after %d iterations\n\n', iter - 1);
    fprintf('%-6s %14s %14s\n', 'Pipe', 'Flow (L/s)', 'Head loss (m)');
    for i = 1:9
        fprintf('Q%-5d %14.2f %14.4f\n', i, Q(i) * 1000, h_loss(i));
    end
    fprintf('\n%-6s %18s\n', 'Unit', 'Pressure loss (kPa)');
    for i = 1:6
        fprintf('S%-5d %18.4f\n', i, p_loss(i));
    end
end

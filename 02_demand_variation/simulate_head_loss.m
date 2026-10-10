function [Q, head_loss] = simulate_head_loss(S1, S2, S3, S4, S5, S6)
% SIMULATE_HEAD_LOSS  Network solve for one set of unit demands S1..S6.
%   Same network as problem 1, written in terms of pipe velocity: the
%   friction factor is updated from the Reynolds number each iteration
%   (laminar, f = 64/Re) until the flows converge.
    Length   = [300, 50, 50, 300, 50, 50, 50, 50, 50];      % m
    Diameter = [0.3, 0.1, 0.1, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1]; % m
    Density  = 1000;     % kg/m^3
    Gravity  = 9.81;     % m/s^2
    U = ones(1, 9);      % starting guess for velocity in each pipe (m/s)

    tol = 6e-6;
    err = 1;

    while err > tol
        Re = (Density .* U .* Diameter) / 0.001002;
        f  = 64 ./ Re;
        k  = (8 .* f .* Length) ./ (pi^2 .* Gravity .* Diameter.^6);
        q  = U .* (pi .* Diameter.^2 ./ 4);

        % Linearise h = k*q^2 about the current flows
        B      = 2 .* k .* q;
        Bprime = -k .* q.^2;

        A = [ 1    -1     0    0    1     0     0    0     0;
              0     1     1    0    0    -1     0    0     0;
              0     0    -1    1    0     0     1    0     0;
              0     0     0    0   -1     0     0    1     0;
              0     0     0    0    0     1     0   -1    -1;
              0     0     0    0    0     0    -1    0     1;
            -B(1) -B(2)  B(3) B(4)  0     0     0    0     0;
              0    B(2)   0    0   B(5)  B(6)   0   B(8)   0;
              0     0   -B(3)  0    0   -B(6) -B(7)  0   -B(9)];

        S = [S1; S2; S3; S4; S5; S6;
             Bprime(1) + Bprime(2) - Bprime(3) - Bprime(4);
            -Bprime(2) - Bprime(5) - Bprime(6) - Bprime(8);
             Bprime(3) + Bprime(6) + Bprime(7) + Bprime(9)];

        Q = gaussPiv(A, S)';

        err = max(abs(q - Q));
        U = Q ./ (pi .* Diameter.^2 ./ 4);
    end

    head_loss = (f .* Length ./ Diameter) .* U.^2 / (2 * Gravity);   % Darcy-Weisbach
end

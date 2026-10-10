function x = gaussPiv(A, b)
% GAUSSPIV  Solve A*x = b by Gaussian elimination with partial pivoting.
%   Partial pivoting (swapping in the row with the largest pivot) avoids
%   dividing by small numbers, which keeps round-off error under control.
    n = size(A, 1);
    Ab = [A b];

    % Forward elimination
    for k = 1:n-1
        [~, pivot] = max(abs(Ab(k:n, k)));
        pivot = pivot + k - 1;
        if pivot ~= k
            Ab([k pivot], :) = Ab([pivot k], :);
        end
        for i = k+1:n
            factor = Ab(i, k) / Ab(k, k);
            Ab(i, k:n+1) = Ab(i, k:n+1) - factor * Ab(k, k:n+1);
        end
    end

    % Back substitution
    x = zeros(n, 1);
    x(n) = Ab(n, end) / Ab(n, n);
    for i = n-1:-1:1
        x(i) = (Ab(i, end) - Ab(i, i+1:n) * x(i+1:n)) / Ab(i, i);
    end
end

function x = gaussPiv(A, b)
    % Get the number of equations
    n = size(A, 1);

    % Augmented matrix
    Ab = [A b];

    % Forward elimination with partial pivoting
    for k = 1:n-1
        % Find the pivot row
        [~, pivot] = max(abs(Ab(k:n, k)));
        pivot = pivot + k - 1;
        
        % Swap rows if necessary
        if pivot ~= k
            Ab([k pivot], :) = Ab([pivot k], :);
        end
        
        % Eliminate entries below the pivot
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
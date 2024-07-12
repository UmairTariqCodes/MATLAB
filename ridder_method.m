function root = ridder1(func, x1, x2, tol)
% Ridder's method for computing the root of f(x) = 0
% USAGE: root = ridder(func, a, b, tol)
% INPUT:
% func = handle of function that returns f(x).
% x1, x2  = limits of the interval containing the root.
% tol  = error tolerance (default is 1.0e6*eps).
% OUTPUT:
% root = zero of f(x) (root = NaN if failed to converge).

if nargin < 4; tol = 1.0e6*eps; end

f1 = func(x1);
if f1 == 0; root = x1; return; end

f2 = func(x2);
if f2 == 0; root = x2; return; end

if sign(f1) == sign(f2)
    error('Root is not bracketed in (%f, %f). f(x1) = %f, f(x2) = %f', x1, x2, f1, f2);
end

for i = 0:30
    % Compute improved root from Ridder's formula
    x3 = 0.5 * (x1 + x2); 
    f3 = func(x3);
    if f3 == 0; root = x3; return; end
    s = sqrt(f3^2 - f1 * f2);
    if s == 0; root = NaN; return; end
    dx = (x3 - x1) * f3 / s;
    if (f1 - f2) < 0; dx = -dx; end
    x4 = x3 + dx; 
    f4 = func(x4);
    % Test for convergence
    if i > 0
        if abs(x4 - xOld) < tol * max(abs(x4), 1.0)
            root = x4; return
        end
    end
    xOld = x4;
    % Re-bracket the root
    if sign(f3) == sign(f4)
        if sign(f1) ~= sign(f4); x2 = x4; f2 = f4;
        else          x1 = x4; f1 = f4;
        end
    else
        x1 = x3; x2 = x4; f1 = f3; f2 = f4;
    end
end

root = NaN;

end

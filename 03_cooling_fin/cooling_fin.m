% COOLING_FIN  2D transient temperature field in a cooling fin (explicit finite differences).
%   Forward Euler in time, second-order central differences in space.
%   Ghost nodes outside the domain impose the boundary conditions:
%     left   - constant heat flux in      right  - adiabatic
%     top    - convection to ambient      bottom - adiabatic (symmetry line)
clear; close all;

h   = 938;      % convective heat transfer coefficient (W/m^2K)
q   = 8312.67;  % heat flux into the base (W/m^2)
L   = 0.04;     % fin length (m)
th  = 0.01;     % fin thickness (m)
k   = 1;        % thermal conductivity (W/mK)
rho = 2100;     % density (kg/m^3)
c   = 850;      % specific heat (J/kgK)
TAmb = 20;      % ambient temperature (C)

M  = 100;       % nodes along x
N  = 20;        % nodes across y
dx = L / M;
dy = th / N;
alpha = k / (rho * c);

% Explicit scheme is stable only if alpha*dt*(1/dx^2 + 1/dy^2) <= 1/2
dt_max = 0.5 / (alpha * (1/dx^2 + 1/dy^2));
dt = 0.001;
fprintf('Stability limit dt = %.4f s, using dt = %.4f s\n', dt_max, dt);

T = TAmb * ones(M + 2, N + 2);   % includes one ring of ghost nodes
T_new = T;
num_timesteps = 10000;

for n = 1:num_timesteps
    for i = 2:M+1
        for j = 2:N+1
            T_new(i,j) = T(i,j) + dt * alpha * ...
                ((T(i+1,j) - 2*T(i,j) + T(i-1,j)) / dx^2 + ...
                 (T(i,j+1) - 2*T(i,j) + T(i,j-1)) / dy^2);
        end
    end

    T_new(2:M+1, 1)   = T_new(2:M+1, 2);                                          % bottom: adiabatic
    T_new(2:M+1, N+2) = T_new(2:M+1, N+1) + h * dy / k * (TAmb - T_new(2:M+1, N+1)); % top: convection
    T_new(1, 2:N+1)   = T_new(2, 2:N+1) + q * dx / k;                             % left: heat flux
    T_new(M+2, 2:N+1) = T_new(M+1, 2:N+1);                                        % right: adiabatic

    T = T_new;
    if mod(n, 1000) == 0
        fprintf('Step %d of %d\n', n, num_timesteps);
    end
end
fprintf('Simulated %.1f s. Peak temperature %.2f C\n', num_timesteps * dt, max(max(T(2:end-1, 2:end-1))));

figure;
imagesc(T(2:end-1, 2:end-1)');
axis xy; colorbar;
title(sprintf('Temperature field after %.0f s (C)', num_timesteps * dt));
xlabel('Node along length'); ylabel('Node across thickness');

% Mirror about the symmetry line to show the full fin cross-section
half = rot90(T(2:end-1, 2:end-1));
figure;
contourf([half; flip(half)], 150, 'EdgeColor', 'none');
colorbar;
title('Temperature contours, full fin section (C)');
xlabel('Node along length'); ylabel('Node across thickness');

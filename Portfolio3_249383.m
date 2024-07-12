clc;
close all;
clear all;

% Parameters
% Heat transfer coefficient (W/m^2-K)
h = 938; % W/m^2-K (3rd, 4th, and 5th digits of your candidate number)

% Dimensions of the cooling fin
L = 0.04; % m, length of the fin
th = 0.01; % m, thickness of the fin

% Thermal properties of the material
k = 1; % W/m-K, thermal conductivity
rho = 2100; % kg/m^3, density
c = 850; % J/kg-K, specific heat capacity

% Discretisation parameters
M = 100; % Number of nodes in the x-direction
N = 20; % Number of nodes in the y-direction
dt = 0.001; % Time step (s), adjusted for stability
dx = L / M; % Mesh size in x-direction (m)
dy = th / N; % Mesh size in y-direction (m)

% Thermal diffusivity
alpha = k / (rho * c);

% Ambient temperature (°C)
TAmb = 20;

% Heat flux (W/m^2) based on candidate number
q = 8312.67; % based on candidate number 249383

% Initialisation
% Initialise the temperature matrix with ambient temperature
% The matrix has extra rows and columns to represent ghost nodes
T = TAmb * ones(M + 2, N + 2);

% Temporary matrix to store new temperature values
T_new = T;

% Simulation parameters
num_timesteps = 10000; % Total number of time steps
save_interval = 1000; % Interval for saving and displaying progress

% Main Simulation Loop
for n = 1:num_timesteps
    % Step 1: Update internal nodes
    for i = 2:M+1 % Loop over internal nodes in x-direction
        for j = 2:N+1 % Loop over internal nodes in y-direction
            % Update temperature using finite difference method
            T_new(i,j) = T(i,j) + dt * alpha * ...
                ((T(i+1,j) - 2*T(i,j) + T(i-1,j)) / dx^2 + ...
                (T(i,j+1) - 2*T(i,j) + T(i,j-1)) / dy^2);
        end
    end
    
    % Step 2: Apply boundary conditions
    % Bottom boundary (adiabatic: no heat transfer)
    T_new(2:M+1,1) = T_new(2:M+1,2);
    
    % Top boundary (convective heat transfer)
    T_new(2:M+1,N+2) = T_new(2:M+1,N+1) + h * dy / k * (TAmb - T_new(2:M+1,N+1));
    
    % Left boundary (constant heat flux)
    T_new(1,2:N+1) = T_new(2,2:N+1) + q * dx / k;
    
    % Right boundary (adiabatic: no heat transfer)
    T_new(M+2,2:N+1) = T_new(M+1,2:N+1);

    % Step 3: Update the temperature matrix
    T = T_new;
    
    % Step 4: Save and display progress at specified intervals
    if mod(n, save_interval) == 0
        disp(['Time step: ', num2str(n)]);
        intermediate_filename = sprintf('temperature_step_%d.mat', n);
        save(intermediate_filename, 'T');
    end
end

% Final Display
disp('Simulation complete.');

% Plotting the Final Temperature Distribution
figure;
% Plot temperature distribution excluding ghost nodes
imagesc(T(2:end-1, 2:end-1));
colorbar;
title('Temperature Distribution at Final Time Step');
xlabel('X Position');
ylabel('Y Position');

% Plotting the Final Contour Plot
% Optional mirroring for better visualisation
mirror = flip(rot90(T(2:end-1, 2:end-1)));
mirror = [rot90(T(2:end-1, 2:end-1)); mirror];

figure;
hold on;
% Create contour plot of the mirrored temperature data
contourf(mirror, 150, 'EdgeColor', 'none');
colorbar;
title('Temperature Contour Plot');
xlabel('X Position');
ylabel('Y Position');


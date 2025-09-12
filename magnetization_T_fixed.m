%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo Simulation of the 2D Ising Model (Metropolis algorithm)
%
% Author: Giuseppe Francesco Conte
% Date:   [29/02/2025]
% Version: 1.0
%
% Description:
%   - Coupling constant J = 1
%   - Lattice size: N x N
%   - Periodic boundary conditions
%   - Temperature T fixed
%   - Random initial configuration with spins in {−1, +1}
%   - Metropolis rule for spin updates
%   - Magnetization computed after each sweep
%   - Snapshots shown at selected sweeps and at the end
%
% Notes:
%   - T_c = 2 / log(1 + sqrt(2)) (critical temperature in 2D Ising model)
%   - The script plots magnetization vs sweeps and a zoomed inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

J = 1; % Coupling constant
N = 100; % Number of spins
Nflip_tentati = 1000000; % Number of attempted spin flips (total sweeps)
T_c = 2/log(1+sqrt(2))
T = 1.5 ; % Fixed temperature
% Generate a random N×N spin configuration in {−1, +1}
M = rand(N);
M(find(M < 0.5)) = -1;
M(find(M > 0.5)) = 1;

% Display the initial random configuration
imagesc(M)
% Custom colormap
colormap([0 0 1; 1 1 0]); % Blue (cold): -1, Yellow (hot): +1
% Set color limits
caxis([-1 1])
% Add colorbar
cb = colorbar;
cb.Ticks = [-1, 1]; % Set colorbar ticks
cb.TickLabels = {'-1', '+1'}; % Custom tick labels
% Image title
title(['Snapshot for $N_{sweeps}$ = 0 '],'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
figure;

% Starting arrays ----> preallocation
M_t = zeros(1, Nflip_tentati);

for ns = 1:Nflip_tentati % Attempted sweeps loop
    
    for n = 1:N^2 % Flip N^2 spins per sweep
    
        % Select a random spin on the grid
        i = randi(N);
        j = randi(N);
    
        % Proposed flip for the selected spin
        M_flip = -M(i,j);
    
        % Periodic boundary neighbors
        if i == 1
            left = N;
            right = i + 1;
        elseif i == N
            left = i - 1;
            right = 1;
        else
            left = i - 1;
            right = i + 1;
        end
    
        if j == 1
            up = N;
            down = j + 1;
        elseif j == N
            up = j - 1;
            down = 1;
        else
            up = j - 1;
            down = j + 1;
        end
    
        % Energy change for the proposed flip
        deltaH = - J * (M_flip - M(i,j)) * (M(i,down) + M(i,up) + M(left,j) + M(right,j));
    
        % Metropolis acceptance criterion
        if deltaH < 0 % If ΔH < 0, accept the flip
            M(i,j) = M_flip;
        elseif deltaH == 0 % If ΔH = 0, accept with probability 1/2
            if rand < 0.5
                M(i,j) = M_flip;
            end
        else % If ΔH > 0, accept with probability exp(-ΔH/T)
            if rand < exp(-deltaH/T)
                M(i,j) = M_flip;
            end
        end
    
    end
    
    % Show intermediate snapshots
    if ns == 10
        imagesc(M)
        % Custom colormap (same as above)
        colormap([0 0 1; 1 1 0]); % Blue: -1, Yellow: +1
        % Set color limits
        caxis([-1 1])
        % Colorbar with ticks at -1 and +1
        cb = colorbar;
        cb.Ticks = [-1, 1]; % Set colorbar ticks
        cb.TickLabels = {'-1', '+1'}; % Custom tick labels
        % Add image title
        title(['Snapshot for $N_{sweeps}$ = ', num2str(ns)],'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        figure;
    end    



    if ns == 100
        imagesc(M)
        % Custom colormap (same as above)
        colormap([0 0 1; 1 1 0]); % Blue: -1, Yellow: +1
        % Set color limits
        caxis([-1 1])
        % Colorbar with ticks at -1 and +1
        cb = colorbar;
        cb.Ticks = [-1, 1]; % Set colorbar ticks
        cb.TickLabels = {'-1', '+1'}; % Custom tick labels
        % Add image title
        title(['Snapshot for $N_{sweeps}$ = ', num2str(ns)],'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
        figure;
    end  


    % Evaluate magnetization at fixed temperature AFTER completing the sweep
    M_m = sum(sum(M)) / N^2;% Store average magnetization of the configuration
    M_t(ns) = abs(M_m);  % Store absolute average magnetization of the configuration
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    AFTER COMPLETING ALL THE SWEEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imagesc(M)
% Final snapshot with the same custom colormap
colormap([0 0 1; 1 1 0]); % Blue: -1, Yellow: +1
% Set color limits
caxis([-1 1])
% Colorbar with ticks at -1 and +1
cb = colorbar;
cb.Ticks = [-1, 1]; % Set colorbar ticks
cb.TickLabels = {'-1', '+1'}; % Custom tick labels
% Add image title
title(['Snapshot for $N_{sweeps}$ = ', num2str(Nflip_tentati)],'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
% Plot M_t vs N_sweeps
plot(1:Nflip_tentati, M_t, '-', 'LineWidth', 2);
ylim([0,0.4])
xlabel('$N_{sweep}$', 'Interpreter', 'latex');
ylabel('$|M|$', 'Interpreter', 'latex');
title(['Absolute Magnetization vs Sweeps at T = ', num2str(T)], 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

 % Add an inset plot
axes('Position', [0.45 0.45 0.4 0.4]); % Inset position [x y width height]
box on; % Add a border around the inset
plot(1:300000, M_t(1:300000), '-', 'LineWidth', 1.5); % Plot the first 300000 steps (note: title below says 100000)
xlabel('$N_{sweep}$', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$|M|$', 'Interpreter', 'latex', 'FontSize', 10);
title('First 100000 Sweeps', 'FontSize', 10); % Title text does not match the 300000 range
grid on;

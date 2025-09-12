%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo simulation of the 2D Ising model (Metropolis algorithm)
% Magnetization as a function of temperature
%
% Difference with the other code:
% Here the temperature is scanned over a range and the average equilibrium
% magnetization is computed and compared with Onsager's exact solution.
% The other code keeps T fixed and shows the time evolution of |M| with sweeps.
%
% Author: [Your Name]
% Date:   [DD/MM/YYYY]
% Version: 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

J = 1; % Coupling constant
N = 50; % Number of spins
Nconf = 60; % Total number of configurations to generate (temperature points)
Nsweeps = 1000000; % Number of sweeps for each temperature value (large value reduces noise at low T)
tau = Nsweeps/5 % Thermalization time
conf = {};
M_average_time=0;
T_c = 2/log(1+sqrt(2)) % Critical temperature in 2D Ising model
temp = linspace(0.1, 4, Nconf); % Temperature range

% Generate an initial random configuration (+1 / -1)
M = rand(N);
M(find(M < 0.5)) = -1;
M(find(M > 0.5)) = 1;
M_time = zeros(1, Nsweeps);
M_t = zeros(1, Nconf);  % Preallocate array for magnetization

for nc = 1:Nconf % Loop over temperatures
    T = temp(nc); % Current temperature

    % Perform Nsweeps sweeps at this temperature
    for ns = 1:Nsweeps % Loop over sweeps
        
        for n = 1:N^2 % Attempt N^2 spin flips per sweep
    
            % Select a random spin
            i = randi(N);
            j = randi(N);
    
            % Flip proposal for the chosen spin
            M_flip = -M(i,j);
    
            % Periodic boundary conditions
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
    
            % Energy difference after the flip
            deltaH = - J * (M_flip - M(i,j)) * (M(i,down) + M(i,up) + M(left,j) + M(right,j));
    
            % Metropolis rule
            if deltaH < 0 % Accept if ΔH < 0
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
         % Magnetization per spin after the sweep
         M_m = sum(sum(M)) / N^2;  
         M_time(ns)= abs(M_m);
    end
    
    % Average magnetization after equilibration
    M_average_time = mean(M_time(tau:end));
    M_t(nc) = abs(M_average_time);  
    
    % Store configuration for this temperature
    conf{nc} = reshape(M, N^2, 1); 
end

% Theoretical magnetization (Onsager solution)
M_theory = zeros(size(temp));
for i = 1:length(temp)
    if temp(i) < T_c
        M_theory(i) = (1 - sinh(2/temp(i))^(-4))^(1/8);
    else
        M_theory(i) = 0;  % No magnetization above T_c
    end
end

% Plot simulated vs theoretical magnetization
figure;
plot(temp, M_t, 'o', 'LineWidth', 2, 'DisplayName', 'Simulated Magnetization');
hold on;
plot(temp, M_theory, '-', 'LineWidth', 2, 'DisplayName', 'Theoretical Magnetization');
xlabel('Temperature (T)');
ylabel('Magnetization (|M|)');
title('Magnetization vs Temperature');
legend('Location', 'Best');  
grid on;
hold off;

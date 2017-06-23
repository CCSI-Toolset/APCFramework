u_plant = [2.46; 30.9; 3.0; 30.9]; % -- [Q1; Q3; Q5; Q7]

% y_plant = [6.9251; 6.9251]; % -- [pH1; pH2]

fprintf('\n');
fprintf('Please wait while the plant is invoked and run for some time-period to extract \n');
fprintf('nominal values corresponding to following input values (u_plant)... \n');
disp(u_plant);
Ts_plant = 0.1; % -- 20 x 0.005 (20 time-steps)
t_plant = 0;
SimulatePlant;
fprintf('-- Nominal output values (y_plant) generated -- \n');
disp(y_plant);

u_name = {'Q_1'; 'Q_3'; 'Q_5'; 'Q_7'};
y_name = {'pH_1'; 'pH_2'};
u_unit = {'m^3/hr'; 'm^3/hr'; 'm^3/hr'; 'm^3/hr'};
y_unit = {''; ''};
t_unit = 'hr';

u_plant = [600000; 6750]; % -- [F_sorbent; F_fluegas]

%y_plant = [87.7645; 600000; 6750]; % -- [CO2capture; F_sorbent_PV; F_fluegas_PV]

fprintf('\n');
fprintf('Please wait while the plant is invoked and run for some time-period to extract \n');
fprintf('nominal values corresponding to following input values (u_plant)... \n');
disp(u_plant);
Ts_plant = 200; % -- 20 x 10 (20 time-steps)
t_plant = 0;
SimulatePlant;
fprintf('-- Nominal output values (y_plant) generated -- \n');
disp(y_plant);

u_name = {'F_{sorbent}'; 'F_{fluegas}'};
y_name = {'CO_2 Capture'; 'F_{sorbent} PV'; 'F_{fluegas} PV'};
u_unit = {'kg/hr'; 'kmol/hr'};
y_unit = {'%'; 'kg/hr'; 'kmol/hr'};
t_unit = 'seconds';

% u_name = {'F_s'; 'F_f'};
% y_name = {'CO_2 Capture'; 'F_s PV'; 'F_f PV'};
% u_unit = {'kg/hr'; 'kmol/hr'};
% y_unit = {'%'; 'kg/hr'; 'kmol/hr'};
% t_unit = 'seconds';

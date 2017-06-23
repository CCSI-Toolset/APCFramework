pH_Neut_Plant = pH_Reactor;
u_plant = [2.46; 30.9; 3.0; 30.9]; % -- [Q1;Q3;Q5;Q7]
pH_Neut_Plant.initialize(u_plant);
y_plant = pH_Neut_Plant.y;

u_name = pH_Neut_Plant.u_name;
y_name = pH_Neut_Plant.y_name;
u_unit = pH_Neut_Plant.u_unit;
y_unit = pH_Neut_Plant.y_unit;
t_unit = pH_Neut_Plant.t_unit;

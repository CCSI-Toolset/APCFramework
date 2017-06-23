% CCSI APC Framework
% APC Framework Simulation Tool

% Code written and developed by:
% Priyadarshi Mahapatra, PhD
% URS Corporation / National Energy Technology Laboratory
% 01 June 2015

clc;
rng_seed = 17;
if  verLessThan('matlab','7.13') % Pre R2011a when 'rng' didn't exist
    rng_legacy(rng_seed);
else
    rng(rng_seed);
end
clear rng_seed;

fprintf('\n');
fprintf(' APC Simulation Tool v%1.1f \n', APCFrameworkVersion('simulation'));
fprintf('-------------------------- \n');
fprintf('This tool provides a SIMULATION platform for the APC Framework. \n');
fprintf('The high-fidelity simulation (''plant'') and APC (''controller'') object \n');
fprintf('need to be setup prior to using this tool. Please exit out of the session \n');
fprintf('(Ctrl-C), in case the plant ''marching'' script or APC object(s) are not \n');
fprintf('defined earlier. Please refer to documentation, with specific examples, for \n');
fprintf('information on how to configure them. \n');
fprintf('\n');
fprintf('Press any key to continue...'); pause;
fprintf('\n\n');

% %% Provide 'Plant' Information
% % Comment - Using run(SimPlantFile) slows down the simulation. Avoid using it.
% 
% fprintf('1. Specify the high-fidelity (''plant'') simulation script: \n');
% plant_specified = 0;
% while ~plant_specified
%     fprintf('   Press any key to browse the script file...'); pause;
%     [FileName,PathName] = uigetfile('*.m','Select the plant-simulation script file');
%     fprintf('\n');
%     if FileName == 0
%         fprintf('     Error: You must provide a script file to continue or press Ctrl-C to exit.\n');
%     else
%         plant_specified = 1;
%     end
% end
% SimPlantFile = [PathName FileName];
% fprintf('-- Plant script ''%s'' specified --\n',SimPlantFile);
% clear FileName PathName plant_specified;
% fprintf('\n');

if ~exist('u_plant','var') || ~exist('y_plant','var')
    fprintf('ERROR: High-fidelity (''plant'') has NOT been previously instantiated. Variables ''u_plant'' \n');
    fprintf('  and ''y_plant'' must be present in the workspace. Please see documentation for guidelines \n');
    fprintf('  on how to define/initialize a ''plant''. Press any key to exit...'); pause;
    fprintf('\n\n');
    return;
end
nup = length(u_plant); % No. of input(s) in high fidelity model
nyp = length(y_plant); % No. of output(s) in high fidelity model

%% Provide Controller Information
fprintf('1. Controller (APC Framework object) specification: \n');
fprintf('   The variable named ''APC'' must be loaded into workspace with the relevant \n');
fprintf('   control-configuration. ');
in = input('Is the correct APC object already present in workspace? (y/n): ','s');
if ~(strcmpi(in,'y'))
    APC_unspecified = 1;
    while APC_unspecified
        fprintf('   Press any key to browse the MAT file containing APC object...'); pause;
        [FileName,PathName] = uigetfile('*.mat','Select the data file containing APC object');
        fprintf('\n');
        if isequal(FileName,0)
            fprintf('     Error: You must provide a MAT file to continue or press Ctrl-C to exit. \n');
        else
            APC_unspecified = 0;
        end
    end
    load([PathName FileName]);
    fprintf('-- Loaded file into workspace -- \n');
end
fprintf('-- Reading %s controller from workspace -- \n',class(APC));
clear FileName PathName APC_unspecified in;
fprintf('\n');

% Some nomenclature
if isa(APC,'APC_MMPC')
    ModL = APC.ModL;
    Mod_idx = APC.Mod_idx;
else
    ModL = 1;
    Mod_idx = 1;
end

Ts = APC.DRM(Mod_idx(1)).Ts;
fprintf('-- Sampling-time is autodetected as %s -- \n',num2str(Ts));
fprintf('\n');

Tend = input('2. Enter the total SIMULATION TIME: ');
fprintf('\n');

fprintf('3. Enter the agitation time-interval. Typically, this is the plant''s process-time / settling time. \n');
fprintf('   For pH Reactor, this is ~0.25 [hr] and for BFB it is ~300 [s]. \n ');
dt_agitation = input('  In what interval should be plant be agitated : ');
N_agitation_gap = round(dt_agitation/Ts);
clear dt_agitation
fprintf('\n');

r_sd = input('4. Enter the standard deviation (in %) of set-point agitation(s): ');
r_sd = r_sd/100;
fprintf('\n');
l_sd = input('5. Enter the standard deviation (in %) of disturbance agitation(s): ');
l_sd = l_sd/100;
fprintf('\n');

fprintf('-- The following provides the CONTROLLER I/O PORTS among plant''s physical ports -- \n');
fprintf('   (Please note the order in which they appear...) \n');
fprintf(['   Manipulated Inputs (MV) - [' num2str(reshape(APC.u_control_idx,1,length(APC.u_control_idx))) '] \n']);
fprintf(['   Controlled Outputs (CV) - [' num2str(reshape(APC.y_control_idx,1,length(APC.y_control_idx))) '] \n']);
fprintf('\n');

l_unspecified = 1;
fprintf('6. Which disturbance(s) you wish to agitate during the simulation, among all plant inputs? \n');
fprintf('   Note - If this interferes with manipulated input selection, disturbance(s) will be imposed \n');
fprintf('          on the corresponding input(s) before control-moves are calculated. Provide empty \n');
fprintf('          vector ''[]'' for no measured disturbance(s). \n');

while l_unspecified
    l_idx = input('   Enter the disturbance(s) index in array format (e.g. [1 3]): ');
    if ~all(ismember(l_idx,1:nup))
        fprintf(['     Error: Indices should be within plant''s input(s) - [' num2str(1:nup) ']. \n']);
    else
        l_unspecified = 0;
    end
end
clear l_unspecified
fprintf('\n');
fprintf('Press any key to start the simulations...'); pause;
fprintf('\n\n');

%% Initialize controller
tic;
APC.initializeController(u_plant);
TimerInitController = toc;

%% Create APC Historian Database
try
    APCResponse = APCHistorian(APC,y_name,u_name,y_unit,u_unit,t_unit);
catch
    APCResponse = APCHistorian(APC);
end
APCResponse.showDynPlot('FrameDuration',Tend/3,'FixedFrameInterval',true,'ShowIntegralError',true,'ShowControlCalcTime',true);

if isa(APC,'APC_MMPC')
    ym = cell(1,ModL);
    for jj = 1:ModL
        ii = Mod_idx(jj);
        ym{jj} = APC.DRM(ii).y;
    end
    APCResponse.saveResponse('Output',y_plant,'Model',ym,'Weight',APC.W_model);
else
    APCResponse.saveResponse('Output',y_plant,'Model',APC.DRM.y,'FilterSD',APC.DRM.SDk);
end

yc_nominal = y_plant(APC.y_control_idx);
r = yc_nominal;
l_nominal = u_plant(l_idx);

%% Run Closed-Loop Simulation
Tbeg = 0;
N = (Tend-Tbeg)/Ts;
t = Tbeg;
TimerTotalSimID = tic;
for k = 1:N
    if rem(k-1-5-N_agitation_gap, 2*N_agitation_gap)==0
        r = yc_nominal.*(1 + r_sd*randn(size(yc_nominal)));
    end
    if rem(k-1-5, 2*N_agitation_gap)==0
        l = l_nominal.*(1 + l_sd*randn(size(l_nominal)));
        u_plant(l_idx) = l;
    end
    
    % Evaluate Controller Move
    tic;
    uc = APC.evalControlMove(r,u_plant,y_plant);
    APCResponse.TimerDB.Control_Calculation(k) = toc;
    APCResponse.TimerDB.Optim_Iteration(k) = APC.countOptimIter;

    % Plant Simulation
    u_plant(APC.u_control_idx) = uc; t_plant = t; Ts_plant = Ts;
    tic;
    SimulatePlant;
    APCResponse.TimerDB.Plant_Simulation(k) = toc;
    
    % Show some friendly progress
    n = 20;
    if rem(k,round(N/n)) == 0
        %fprintf('-- %.0f %% complete --\n',k/N*100);
        fprintf('   [');
        for i = 1:round(k/N*n)
            fprintf('-');
        end
        for i = 1:n-round(k/N*n)
            fprintf(' ');
        end
        fprintf(']\n');
    end
    
    % Save to Historian
    APCResponse.saveResponse('Time',t,'SetPoint',r,'Input',u_plant);
    if isa(APC,'APC_MMPC')
        for jj = 1:ModL
            ii = Mod_idx(jj);
            ym{jj} = APC.DRM(ii).y;
        end
        APCResponse.saveResponse('Output',y_plant,'Model',ym,'Weight',APC.W_model);
    else
        APCResponse.saveResponse('Output',y_plant,'Model',APC.DRM.y,'FilterSD',APC.DRM.SDk);
    end
    
    t = t + Ts;
end
% Save to Historian (R(:,k+1) = R(:,k); U(:,k+1) = U(:,k);)
APCResponse.saveResponse('Time',t,'SetPoint',r,'Input',u_plant);
% Save Timers
APCResponse.TimerDB.Total_Simulation = toc(TimerTotalSimID);
APCResponse.TimerDB.Initialize_Controller = TimerInitController;
% Show Brief Summary on Screen
fprintf('\nTotal Simulation Time = %f min\n',APCResponse.TimerDB.Total_Simulation/60);
fprintf('\nTotal Cumulative Residual = %f\n',APCResponse.evalControllerError());
% Clear Variables
clear Ts Tbeg Tend yc_nominal l_nominal t r l ym uc yc jj ii k i l_idx r_sd l_sd N_agitation_gap N n
clear TimerInitController TimerTotalSimID;

%% Plotting the responses
APCResponse.displayPlot('ControlVars','AllVars','CompTime','OptimIter','ModWeight');
% APCResponse.plotControllerVariables;
% APCResponse.plotAllVariables;
% APCResponse.plotComputationTime;
% APCResponse.plotOptimIteration
% APCResponse.plotModelWeights;

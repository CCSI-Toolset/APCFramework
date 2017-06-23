% CCSI APC Framework
% Script for generating Step-Response Data (Training Data) using random
% input sequence of specified standard-deviation every fixed specified 
% time-step

% Code written and developed by:
% Priyadarshi Mahapatra
% URS Corporation / National Energy Technology Laboratory / U.S. Department of Energy
% Morgantown, WV 26505
% 20 Nov 2013
%
% -------------------------------------------------------------------------
% Specific instructions for pH_Reactor example -
%
% Make sure pH_Neut_Plant is defined and initialized. If not, use the
% following procedure OR just run the 'Initialize_pH_Reactor.m' script -
% 1. Create a pH_Reactor object:
%        pH_Neut_Plant = pH_Reactor;
% 2. Specify the nominal input(s):
%        u_plant = [2.46; 30.9; 3.0; 30.9]; % -- [Q1;Q3;Q5;Q7]
% 3. Initialize the reactor object:
%        pH_Neut_Plant.initialize(u_plant); 
% 4. Obtain the nominal outputs(s):
%        y_plant = pH_Neut_Plant.y;
% 5. Specify the I/O names and units (optional):
%        u_name = pH_Neut_Plant.u_name;
%        y_name = pH_Neut_Plant.y_name;
%        u_unit = pH_Neut_Plant.u_unit;
%        y_unit = pH_Neut_Plant.y_unit;
%        t_unit = pH_Neut_Plant.t_unit;

%% Specification set for pH-Reactor example 
% If you want to skip the steps in the script, uncomment the section below and 
% copy-paste in command window

% pH_Neut_Plant = pH_Reactor;
% u_plant = [2.46; 30.9; 3.0; 30.9]; % -- [Q1;Q3;Q5;Q7]
% pH_Neut_Plant.initialize(u_plant);
% y_plant = pH_Neut_Plant.y;
% ModL = 3;
% Tend{1} = 10;
% Tend{2} = 10;
% Tend{3} = 10;
% Ts{1} = 0.005;
% Ts{2} = 0.005;
% Ts{3} = 0.005;
% um_idx{1} = [2 4 1]; ym_idx{1} = [1 2];
% um_idx{2} = [2 4 1]; ym_idx{2} = [1 2];
% um_idx{3} = [2 4 1]; ym_idx{3} = [1 2];
% u_nominal{1} = [2.46; 30.9; 3.0; 30.9];
% u_nominal{2} = [1.5; 24; 3; 30];
% u_nominal{3} = [4; 40; 3; 25];
% Tend{1} = 10;
% Tend{2} = 10;
% Tend{3} = 10;

%% Display introduction text
clc;
fprintf('\n');
fprintf(' Step-Response Data Generator  \n');
fprintf('------------------------------  \n');
fprintf('This tool helps the user to generate Step-Response Time-Series Data (Training Data) \n');
fprintf('required for developing dynamic reduced order (DRM) model(s). The high-fidelity \n');
fprintf('simulation (''plant'') needs to be defined prior to using this tool. If not, please \n');
fprintf('exit out of the session by pressing Ctrl-C. Please refer to the documentation, with \n');
fprintf('specific examples, for information on how to instantiate a plant. \n');
fprintf('\n');
fprintf('Press any key to continue...'); pause;
fprintf('\n\n');
fprintf('-- Plant is currently initialized as follows --\n');
fprintf('   (transpose ('') denotes column-vector) \n');

if size(u_plant,1) > size(u_plant,2) % denotes column-vector
    fprintf(['   Input(s): [' num2str(u_plant') ']'' \n']);
else
    fprintf(['   Input(s): [' num2str(u_plant) '] \n']);
end
if size(y_plant,1) > size(y_plant,2) % denotes column-vector
    fprintf(['   Output(s): [' num2str(y_plant') ']'' \n']);
else
    fprintf(['   Output(s): [' num2str(y_plant) '] \n']);
end
fprintf('\n');
%% Define number of training-sets
ModL = input('1. Enter the desired number of training-sets to be generated: ');
fprintf('\n');

%% Define sampling time
% Note - This has to be same for all DRM(s)
fprintf('2. Define the sampling time for each training set. \n');
for jj = 1:ModL
    Ts{jj} = input(['   Enter sampling time (in same time-unit as ''plant'') for Set-' num2str(jj) ': ']);
end
clear jj
fprintf('\n');

%% Define I/O variables for training-sets
fprintf('3. Define input/output (I/O) variables among plant I/O variables for each training-set \n');
fprintf('   through variable-indexing. For example, if you wish to have DRM contain inputs [Q3 Q7 Q1] \n');
fprintf('   among plant''s inputs [Q1 Q3 Q5 Q7], the DRM input-index is [2 4 1]. Order of inputs will \n');
fprintf('   be preserved. Please refer to documentation for use-case examples. \n');

for jj = 1:ModL
    um_unspecified = 1;
    while um_unspecified
        um_idx{jj} = input(['   Enter the desired input(s) in array format (e.g. [2 4 1]) for Set-' num2str(jj) ': ']);
        if ~all(ismember(um_idx{jj},1:length(u_plant)))
            fprintf(['     Error: Indices should be within plant''s input(s) - [' num2str(1:length(u_plant)) ']. \n']);
        else
            um_unspecified = 0;
        end
    end
    ym_unspecified = 1;
    while ym_unspecified
        ym_idx{jj} = input(['   Enter the desired output(s) in array format (e.g. [1 2]) for Set-' num2str(jj) ': ']);
        if ~all(ismember(ym_idx{jj},1:length(y_plant)))
            fprintf(['     Error: Indices should be within plant''s output(s) - [' num2str(1:length(y_plant)) ']. \n']);
        else
            ym_unspecified = 0;
        end
    end       
end
clear um_unspecified ym_unspecified jj
fprintf('\n');

%% Define the plant's nominal operating point for each training-set
fprintf('4. Define plant''s nominal operating point for each perturbation. Typically choose points \n');
fprintf('   which may span a range of nonlinear operations frequently encountered during load changes \n');
fprintf('   or disturbance rejection scenarios. \n');

for jj = 1:ModL
    up_unspecified = 1;
    while up_unspecified
        u_nominal{jj} = input(['   Enter the nominal plant input(s) in array format (e.g. [4 40 3 25]) for Set-' num2str(jj) ': ']);
        if ~isequal(size(u_nominal{jj}),size(u_plant))
            if size(u_plant,1) > size(u_plant,2) % denotes column-vector
                space_ch = ';';
            else
                space_ch = ',';
            end
            fprintf('     Error: Input array size/orientation does not match plant''s input definition. Remember to \n');
            fprintf(['            separate the input(s) using ''' space_ch ''' character. \n']);
            clear space_ch
        else
            up_unspecified = 0;
        end
    end
end
clear up_unspecified jj
fprintf('\n');

%% Specify number of samples for each training-set
fprintf('5. Define amount of training required for each set by specifying end-time. \n');
for jj = 1:ModL
    fprintf(['   Sampling-time for Set-' num2str(jj) ' was provided earlier as ' num2str(Ts{jj}) ' \n']);
    Tend{jj} = input('   Enter end-time (in same time-unit as ''plant'') for this training set : ');
end
clear jj
fprintf('\n');

%% ======================= USER-INPUTS COMPLETED ==========================
% Run Simulation
% T, U and Y contain the step-response data for each DRM

y_nominal = cell(1,ModL);
T = cell(1,ModL);
U = cell(1,ModL); Y = cell(1,ModL);
Uss = cell(1,ModL); Yss = cell(1,ModL);

for jj = 1:ModL

    T{jj} = 0:Ts{jj}:Tend{jj};

    % Find the nominal output values for each DRM.
    % This is very specific for the pH Reactor example.
    % In many cases such as ACM-based plant, initializing the plant is not trivial.
    pH_Neut_Plant.initialize(u_nominal{jj});
    y_nominal{jj} = y_plant;
    fprintf(['-- Plant initialized at nominal operating point ' num2str(jj) ' --\n']);
    
    Uss{jj} = u_nominal{jj}(um_idx{jj});
    Yss{jj} = y_nominal{jj}(ym_idx{jj});

    % Start with the nominal values
    U{jj}(:,1) = Uss{jj};
    Y{jj}(:,1) = Yss{jj};

    % Run plant step-response simulations
    fprintf(['-- Generating step-responses for DRM-' num2str(jj) '...']);
    for k = 1:length(T{jj})-1
        if k > 1
            U{jj}(:,k) = U{jj}(:,k-1);
            if rem(k-1-10,50)==0 % Trigger change in input u after every 20 time-steps
                % Step-change with standard-deviation of 2% in each input
                %u_plant = u_nominal{jj}.*(1+0.02*randn(size(u_nominal{jj})));
                U{jj}(:,k) = Uss{jj}.*(1+0.02*randn(size(Uss{jj}))); 
            end
        end
        % Plant simulation
        u_plant = u_nominal{jj}; u_plant(um_idx{jj}) = U{jj}(:,k);
        t_plant = T{jj}(k);
        Ts_plant = Ts{jj};
        SimulatePlant;
        Y{jj}(:,k+1) = y_plant(ym_idx{jj});
    end
    U{jj}(:,k+1) = U{jj}(:,k);
    fprintf('done.\n\n');
    
    if exist('u_name','var') &&  length(um_idx{jj}) <= length(u_name)
        um_name{jj} = u_name(um_idx{jj});
    end
    if exist('y_name','var') &&  length(ym_idx{jj}) <= length(y_name)
        ym_name{jj} = y_name(ym_idx{jj});
    end
    if exist('u_unit','var') &&  length(um_idx{jj}) <= length(u_unit)
        um_unit{jj} = u_unit(um_idx{jj});
    end
    if exist('y_unit','var') &&  length(ym_idx{jj}) <= length(y_unit)
        ym_unit{jj} = y_unit(ym_idx{jj});
    end
    if exist('t_unit','var')
        tm_unit{jj} = t_unit;
    end
    
end

clearvars -except U Y T Uss Yss Ts um_idx ym_idx um_name ym_name um_unit ym_unit tm_unit
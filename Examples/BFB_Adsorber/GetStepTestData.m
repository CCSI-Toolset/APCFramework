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
% Specific instructions for BFB_Adsorber example -
%
% Make sure BFBAdsorber simulink model is defined and initialized. If not, 
% the following steps may be used as a guide -
% 1. Setup the ACM (.acmf) file. Make sure the simulation has converged to
%    a steady-state. Take note of the relevant I/O variable names and their
%    nominal values.
% 2. Setup the ACM-embedded Simulink model (.mdl). Associate the ACM-block
%    with the simulation model. Choose the 'plant' input/output variables
%    as noted earlier. Check/rectify number of inputs/outputs in Demux/Mux
%    blocks and make/remove connections as required. Set the communication
%    time (typically same as Ts_plant). Save the model as BFBAdsorber.mdl.

%% Display introduction text
clc;
fprintf('\n');
fprintf(' Step-Response Data Generator  \n');
fprintf('------------------------------  \n');
fprintf('This tool helps the user to generate Step-Response Time-Series Data (Training Data) \n');
fprintf('required for developing dynamic reduced order (DRM) model(s). \n');
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

%% Define training-set no.
jj = input('1. Enter the training-set no. to be generated: ');
fprintf('\n');

%% Define sampling time
% Note - This has to be same for all DRM(s)
Ts{jj} = input('2. Define the sampling time for training (in same time-unit as ''plant''): ');
fprintf('\n');

%% Define I/O variables for training-sets
fprintf('3. Define input/output (I/O) variables among plant I/O variables for each training-set \n');
fprintf('   through variable-indexing. For example, if you wish to have DRM contain inputs [Q3 Q7 Q1] \n');
fprintf('   among plant''s inputs [Q1 Q3 Q5 Q7], the DRM input-index is [2 4 1]. Order of inputs will \n');
fprintf('   be preserved. Please refer to documentation for use-case examples. \n');

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
clear um_unspecified ym_unspecified
fprintf('\n');

%% Define the plant's nominal operating point for each training-set
u_nominal{jj} = u_plant;
y_nominal{jj} = y_plant;

%% Specify number of samples for each training-set
fprintf('5. Define amount of training required for each set by specifying end-time. \n');
fprintf('   Sampling-time for this training-set was provided earlier as %s. \n',num2str(Ts{jj}));
Tend{jj} = input(['   Enter end-time (in same time-unit as ''plant'') for Set-' num2str(jj) ': ']);
fprintf('\n');

%% ======================= USER-INPUTS COMPLETED ==========================
% Run Simulation
% T, U and Y contain the step-response data for each DRM

T{jj} = 0:Ts{jj}:Tend{jj};

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
        if rem(k-1-5,350)==0 % Trigger change in input u after every 20 time-steps
            % Step-change with standard-deviation of 2% in each input
            %u_plant = u_nominal{jj}.*(1+0.02*randn(size(u_nominal{jj})));
            U{jj}(:,k) = Uss{jj}.*(1+0.05*randn(size(Uss{jj}))); 
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

clearvars -except U Y T Uss Yss Ts um_idx ym_idx um_name ym_name um_unit ym_unit tm_unit
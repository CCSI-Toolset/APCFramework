% CCSI APC Framework
% APC Framework Configuration Tool

% Code written and developed by:
% Priyadarshi Mahapatra, PhD
% URS Corporation / National Energy Technology Laboratory
% 13 Sept 2013

%% Display introduction text (I)
clc;
fprintf('\n');
fprintf(' APC Configuration Tool v%1.1f \n', APCFrameworkVersion('config'));
fprintf('----------------------------- \n');
fprintf('This tool provides a CONFIGURATION platform for the APC Framework. \n');
fprintf('The control-model (DRM), generated from high-fidelity (''plant'') simulations, \n');
fprintf('need to be setup prior to using this tool. If the DRM objects are not \n');
fprintf('defined, please exit out of the session by pressing Ctrl-C. Please refer \n');
fprintf('to the documentation, with specific examples, for information on how to \n');
fprintf('configure them. \n');
fprintf('\n');
fprintf('Press any key to continue...'); pause;
fprintf('\n\n');

%% Provide Control-Model(s) / DRM(s) Information
fprintf('I. Essential Controller Specifications \n');
fprintf('-------------------------------------- \n');
fprintf('The following specifications are REQUIRED to instantiate a controller object. \n');
fprintf('\n');
fprintf('1. Control-model(s) / DRM specification: \n');
fprintf('   The variable named ''DRM'' must be loaded into workspace with the relevant \n');
fprintf('   control-model(s). ');
in = input('Is the correct DRM already present in the workspace? (y/n): ','s');
if ~(strcmpi(in,'y'))
    DRM_unspecified = 1;
    while DRM_unspecified
        fprintf('   Press any key to browse the MAT file containing DRM(s)...'); pause;
        fprintf('\n');
        [FileName,PathName] = uigetfile('*.mat','Select the data file containing DRM(s)');
         if FileName == 0
            fprintf('     ERROR: You must provide a MAT file to continue or press Ctrl-C to exit. \n');
        else
            DRM_unspecified = 0;
        end
    end
    load([PathName FileName]);
    fprintf(['-- Loaded DRM file ''' FileName ''' into workspace -- \n']);
    clear FileName PathName DRM_unspecified
end
clear in

try
    nDRM = length(DRM);
catch me
    fprintf('   ERROR: ''DRM'' variable not present in workspace. Press any key to exit...'); pause;
    fprintf('\n\n');
    return;
end

try
    for jj = 1:nDRM
%         % Quick parse of idx and SS to see if compatible
%         DRM(jj).u_idx; DRM(jj).y_idx;
%         DRM(jj).uss; DRM(jj).yss;
        assert(isa(DRM(jj),'DRM_DABNet') || isa(DRM(jj),'DRM_SS'));
    end
    fprintf('-- Reading %d %s model(s) from workspace -- \n',nDRM,class(DRM(1)));
catch me
    fprintf('   ERROR: ''DRM'' empty or incompatible type. Please see documentation for \n');
    fprintf('     ''DRM'' specification. Press any key to exit...'); pause;
    fprintf('\n\n');
    clear nDRM jj
    return;
end
fprintf('\n');

%% Configure Control-Models (among DRM(s)) for MMPC
if isa(DRM,'DRM_SS') && nDRM > 1
    fprintf('** Choose APC "control-model(s)" from within the DRM-bank. Provide in MATLAB array format. \n');
    fprintf('   For example, [1 3] will consider Models 1 & 3 for control purpose in MMPC formulation. \n');
    fprintf('   [ENTER] takes all model(s) in the bank. \n');
    fprintf('   Press any key to quickly go through the nominal value(s) as a reminder...'); pause;
    fprintf('\n\n');

    for jj = 1:nDRM
        fprintf('   -------- DRM %d -------- \n',jj);
        fprintf('   Inputs  : [');
        if isa(DRM(jj).u_name,'cell') && ~isempty(DRM(jj).u_name{1}) && ~isempty(strtrim(DRM(jj).u_name{1}))
            fprintf(DRM(jj).u_name{1});
        else
            fprintf(['u' num2str(DRM(jj).u_idx(1))]);
        end
        for i = 2:length(DRM(jj).u_idx)
            if isa(DRM(jj).u_name,'cell') && ~isempty(DRM(jj).u_name{i}) && ~isempty(strtrim(DRM(jj).u_name{i}))
                fprintf([' ' DRM(jj).u_name{i}]);
            else
                fprintf([' u' num2str(DRM(jj).u_idx(i))]);
            end
        end
        fprintf('] = [');
        fprintf(num2str(DRM(jj).uss(1)));
        for i = 2:length(DRM(jj).uss)
            fprintf([' ' num2str(DRM(jj).uss(i))]);
        end
        fprintf('] \n');
        
        fprintf('   Outputs : [');
        if isa(DRM(jj).y_name,'cell') && ~isempty(DRM(jj).y_name{1}) && ~isempty(strtrim(DRM(jj).y_name{1}))
            fprintf(DRM(jj).y_name{1});
        else
            fprintf(['y' num2str(DRM(jj).y_idx(1))]);
        end
        for i = 2:length(DRM(jj).y_idx)
            if isa(DRM(jj).y_name,'cell') && ~isempty(DRM(jj).y_name{i}) && ~isempty(strtrim(DRM(jj).y_name{i}))
                fprintf([' ' DRM(jj).y_name{i}]);
            else
                fprintf([' y' num2str(DRM(jj).y_idx(i))]);
            end
        end
        fprintf('] = [');
        fprintf(num2str(DRM(jj).yss(1)));
        for i = 2:length(DRM(jj).yss)
            fprintf([' ' num2str(DRM(jj).yss(i))]);
        end
        fprintf('] \n');
    end
    fprintf('\n');

    Mod_unspecified = 1;
    while Mod_unspecified
        Mod_idx = input('   Enter the desired control-model index in array format: ');
        if ~all(ismember(Mod_idx,1:nDRM))
            fprintf('     ERROR: You must provide indices from within 1 to %d.\n', nDRM);
        else
            Mod_unspecified = 0;
        end
    end
    fprintf('\n');
    if isempty(Mod_idx)
        Mod_idx = 1:nDRM;
    end
else
    Mod_idx = 1; % For DABNet or NARMA
end
%fprintf('\n');
fprintf('-- Control-model(s) selection from DRM-bank complete -- \n');
ModL = length(Mod_idx);
clear nDRM Mod_unspecified i j jj

%% Configure Controller Input(s) and Output(s)
fprintf('\n');
fprintf('2. Specify controller input(s)/output(s): \n');
fprintf('   Press any key to continue...'); pause;
fprintf('\n\n');
fprintf('   The selected DRM''s I/O variable(s) are shown below. The indices correspond to the DRM variables \n');
fprintf('   among ''plant'' variables (provided during step-tests). Please refer to documentation, \n');
fprintf('   with specific examples, for more clarity. \n');
fprintf('\n');
for jj = 1:ModL
    ii = Mod_idx(jj);
    um_idx{jj} = DRM(ii).u_idx;
    ym_idx{jj} = DRM(ii).y_idx;
    fprintf('   -------- DRM %d -------- \n',ii);
%     fprintf(['   Inputs  - [' num2str(um_idx{jj}) '] \n']);
%     fprintf(['   Outputs - [' num2str(ym_idx{jj}) '] \n']);
    
    fprintf('   Inputs  - [');
    fprintf(num2str(DRM(ii).u_idx(1)));
    for i = 2:length(DRM(ii).u_idx)
        fprintf([' ' num2str(DRM(ii).u_idx(i))]);
    end
    fprintf('] : [');
    if isa(DRM(jj).u_name,'cell') && ~isempty(DRM(jj).u_name{1}) && ~isempty(strtrim(DRM(jj).u_name{1}))
        fprintf(DRM(jj).u_name{1});
    else
        fprintf(['u' num2str(DRM(jj).u_idx(1))]);
    end
    for i = 2:length(DRM(jj).u_idx)
        if isa(DRM(jj).u_name,'cell') && ~isempty(DRM(jj).u_name{i}) && ~isempty(strtrim(DRM(jj).u_name{i}))
            fprintf([' ' DRM(jj).u_name{i}]);
        else
            fprintf([' u' num2str(DRM(jj).u_idx(i))]);
        end
    end
    fprintf('] \n');

    fprintf('   Outputs - [');
    fprintf(num2str(DRM(ii).y_idx(1)));
    for i = 2:length(DRM(ii).y_idx)
        fprintf([' ' num2str(DRM(ii).y_idx(i))]);
    end
    fprintf('] : [');
    if isa(DRM(jj).y_name,'cell') && ~isempty(DRM(jj).y_name{1}) && ~isempty(strtrim(DRM(jj).y_name{1}))
        fprintf(DRM(jj).y_name{1});
    else
        fprintf(['y' num2str(DRM(jj).y_idx(1))]);
    end
    for i = 2:length(DRM(jj).y_idx)
        if isa(DRM(jj).y_name,'cell') && ~isempty(DRM(jj).y_name{i}) && ~isempty(strtrim(DRM(jj).y_name{i}))
            fprintf([' ' DRM(jj).y_name{i}]);
        else
            fprintf([' y' num2str(DRM(jj).y_idx(i))]);
        end
    end
    fprintf('] \n');
end
fprintf('\n');

[um_available_idx,~] = findArrayIntersection(um_idx);
[ym_available_idx,~] = findArrayIntersection(ym_idx);

fprintf('   The available control-input(s)/output(s) are given as follows. \n');
fprintf(['     Available Control Inputs  - [' num2str(um_available_idx) '] \n']);
fprintf(['     Available Control Outputs - [' num2str(ym_available_idx) '] \n']);
fprintf('\n');

uc_unspecified = 1;
while uc_unspecified
    uc_specified_idx = input('   Enter the desired manipulated input(s) in array format (e.g. [2 4]): ');
    if ~all(ismember(uc_specified_idx,um_available_idx))
        fprintf('     ERROR: You must provide indices corresponding to the available input(s).\n');
    else
        uc_unspecified = 0;
    end
end
if isempty(uc_specified_idx)
    uc_specified_idx = um_available_idx;
end
yc_unspecified = 1;
while yc_unspecified
    yc_specified_idx = input('   Enter the desired control output(s) in array format (e.g. [1 2]): ');
    if ~all(ismember(yc_specified_idx,ym_available_idx))
        fprintf('     ERROR: You must provide indices corresponding to the available output(s).\n');
    else
        yc_unspecified = 0;
    end
end
if isempty(yc_specified_idx)
    yc_specified_idx = ym_available_idx;
end
for jj = 1:ModL
    [~,uc_idx{jj},~] = intersect(um_idx{jj},uc_specified_idx);
    [~,yc_idx{jj},~] = intersect(ym_idx{jj},yc_specified_idx);
    %lc_specified_idx = um_idx{jj};
    %lc_specified_idx(uc_idx{jj}) = [];
    %[~,lc_idx{jj},~] = intersect(um_idx{jj},lc_specified_idx);
end

fprintf('\n');
fprintf('   Note: Remaining inputs in each DRM not specified as ''manipulated-inputs(s)'' will be \n');
fprintf('         designated as ''measured-disturbance(s)'' in controller formulation. \n'); 
fprintf('\n');
fprintf('-- Controller I/O specification complete -- \n');
fprintf('   The following control variable(s) were provided - \n');
fprintf(['     Manipulated Input(s) - [' num2str(um_idx{1}(uc_idx{1})) '] \n']);
fprintf(['     Controlled Output(s) - [' num2str(ym_idx{1}(yc_idx{1})) '] \n']);
clear um_available_idx ym_available_idx uc_specified_idx yc_specified_idx lc_specified_idx uc_unspecified yc_unspecified jj ii
clear ModL um_idx ym_idx

%% Initialize Controller
clear APC
fprintf('\n');
fprintf('   Based on DRM(s) selected, ');
if isa(DRM,'DRM_SS')
    if length(Mod_idx) > 1
        fprintf('Multiple Model Predictive Controller (MMPC)');
        APC = APC_MMPC(DRM,Mod_idx,yc_idx,uc_idx);
    else
        fprintf('Model Predictive Controller (MPC) [Single Model]');
        APC = APC_MPC(DRM(Mod_idx),yc_idx{1},uc_idx{1});
    end
elseif isa(DRM,'DRM_DABNet')
    fprintf('DABNet-based Nonlinear Model Predictive Controller (DABNet-NMPC)');
    APC = APC_NMPC(DRM(Mod_idx),yc_idx{1},uc_idx{1});
elseif isa(DRM,'DRM_NARMA')
    fprintf('NARMA-based Nonlinear Model Predictive Controller (NARMA-NMPC)');
    APC = APC_NMPC(DRM(Mod_idx),yc_idx{1},uc_idx{1});
end
fprintf(' is auto-detected. \n');
fprintf('\n');
fprintf('-- Controller object ''APC'' defined -- \n')
fprintf('\n');
fprintf('Press any key to continue...'); pause;
fprintf('\n\n');

%% Display introduction text (II)
clc;
fprintf('II. Optional Controller Specifications \n');
fprintf('-------------------------------------- \n');
fprintf('The following specifications are not essential at this point and may be configured later. \n');
fprintf('These include definition of following variables - \n');
fprintf('\n');
fprintf('- Input Absolute and Rate Constraints (umin, umax, dumin, dumax) \n');
fprintf('- Output [Hard/Soft] Constraints (ymin, ymax) \n');
fprintf('- Input Weights (wu) \n');
fprintf('- Output Weights (wy) \n');
fprintf('- Prediction Horizon (P) \n');
fprintf('- Control Horizon (M) \n');
fprintf('- Model Covariance Coefficient (lambda) \n');
fprintf('\n');
%fprintf('In newer versions of MATLAB (post R2007), you may define/edit these values directly \n');
%fprintf('in the Variable Editor by double-clicking the APC object in workspace. \n');
% ---- Note to developer - Include this in later version where user can
% change APC details in Variable Editor and later invoke setupController().
% without any arguments. Check for existance of each variable. If it is still
% not defined, let user know and continue default values. Remove case 3.
in = input('Do you want to continue configuring these variables? (y/n): ','s');
if ~(strcmpi(in,'y'))
    fprintf('\n');
    fprintf('-- User declined to configure variables at this point -- \n');
    fprintf('\n');
    %fprintf('Next-Step: Either define them through Variable Editor as mentioned above or issue commands such as - \n');
    fprintf('Next-Step: Issue commands such as - \n\n');
    fprintf('    APCParameters.P = 50; \n');
    fprintf('    APCParameters.umin = [0;0]; \n');
    fprintf('    APCParameters.wy = [1 0; 0 2]; \n');
    if isa(APC,'APC_MMPC')
        fprintf('    APCParameters.lambda = 0.02; \n');
    end
    fprintf('    ...etc. \n\n');
    fprintf('After doing this, invoke this command to setup the controller - \n');
    fprintf('\n');
    fprintf('    APC.setupController(APCParameters); \n');
    
    fprintf('\n');
    fprintf('OR to use default values of controller parameters, simply use it without any arguments - \n');
    fprintf('\n');
    fprintf('    APC.setupController(); \n');
    fprintf('\n');
    fprintf('Press any key to exit the APC Config Tool...'); pause;
    fprintf('\n\n');
    clearvars -except APC
    return;
end
fprintf('\n');

% *************** Still part of development

%% Configure Constraints
% Enter input constraints

% Enter input rate-constraints

% Enter output constraints

%% Configure Weights
% Enter input weights

% Enter output weights

%% Configure Prediction and Control Horizon

%APC.setupController(P,M,umin,umax,dumin,dumax,ymin,ymax,wu,wy,lambda);

clearvars -except APC
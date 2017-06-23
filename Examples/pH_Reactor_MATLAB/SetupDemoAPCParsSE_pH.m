% 'APC' object should be present in the workspace
clear APCParameters
%% 2-Input 2-Output
APCParameters.Tag = 'NMPC-DE-C';
APCParameters.Description = 'pH 2-Tank Centralized DABNet-NMPC w/ Disturbance Estimation';
% Input Hard Constraints
APCParameters.umin = [0;0];
APCParameters.umax = [60;60];
% Input-Rate Hard Constraints
APCParameters.dumin = [-1;-1];
APCParameters.dumax = [1;1];
% Output Hard Constraints
APCParameters.ymin = [0;0];
APCParameters.ymax = [14;14];
% Input Weights
APCParameters.wu = 0.5*diag([1;1]);
% Outputs Weights
APCParameters.wy = diag([1;1]);
% Prediction and Control Horizon
APCParameters.P = 20;
APCParameters.M = 6;

APC.setupController(APCParameters);
clear APCParameters

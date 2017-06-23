% 'APC' object should be present in the workspace
clear APCParameters
%% 2-Input 2-Output
APCParameters.Tag = 'MMPC-C';
APCParameters.Description = 'pH 2-Tank Centralized MMPC';
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
% Model Covariance Coefficient where, Lambda(:) = lambda*wy
APCParameters.lambda = 0.1;

APC.setupController(APCParameters);
clear APCParameters

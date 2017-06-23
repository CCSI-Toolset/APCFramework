% 'APC' object should be present in the workspace
clear APCParameters
%% 1-Input 1-Output
APCParameters.Tag = 'MMPC-C';
APCParameters.Description = 'BFB Adsorber Centralized MMPC';
% Input Hard Constraints
% APCParameters.umin = 0.5*APC.DRM(APC.Mod_idx(1)).uss(APC.uc_idx{1});
% APCParameters.umax = 1.5*APC.DRM(APC.Mod_idx(1)).uss(APC.uc_idx{1});
APCParameters.umin = [300000];
APCParameters.umax = [900000];
% APCParameters.umin = [10000];
% APCParameters.umax = [1200000];
% Input-Rate Hard Constraints
APCParameters.dumin = 2*[-150000];
APCParameters.dumax = 2*[150000];
% Output Hard Constraints
APCParameters.ymin = [70];
APCParameters.ymax = [95];
% APCParameters.ymin = [-Inf];
% APCParameters.ymax = [Inf];
% Input Weights
APCParameters.wu = 5e-7;
% Outputs Weights
APCParameters.wy = 1e5;
% Prediction and Control Horizon
APCParameters.P = 50;
APCParameters.M = 10;
% Model Covariance Coefficient where, Lambda(:) = lambda*wy
APCParameters.lambda = 1e-4;

APC.setupController(APCParameters);
clear APCParameters

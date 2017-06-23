% CCSI APC Framework
% MATLAB's DABNet-Parameter (from DRM-Builder) to APC Framework's DRM Converter Tool

% Code written and developed by:
% Priyadarshi Mahapatra
% URS Corporation / National Energy Technology Laboratory / U.S. Department of Energy
% Morgantown, WV 26505
% 20 Nov 2013

clear all; clc; close all;

fprintf('\n');
fprintf(' MATLAB''s DABNet-Parameter file to APC Framework compatible DRM converter \n');
fprintf('--------------------------------------------------------------------------- \n');
fprintf('This script converts DRM-Builder''s DABNet-Parameter file into APC Framework''s \n');
fprintf('DRM_DABNet object for use in MPC algorithm. \n');
fprintf('\n');
file_unspecified = 1;
while file_unspecified
    fprintf('Press any key to browse the .m file containing DABNet parameters data...'); pause;
    [FileName,PathName] = uigetfile('*.m','Select the file containing DABNet parameters data');
    fprintf('\n');
    if FileName == 0
        fprintf('  Error: You must provide a .m file to continue or press Ctrl-C to exit.\n');
    else
        file_unspecified = 0;
    end
end
try
    run([PathName FileName]); clear FileName PathName file_unspecified;
catch
    fprintf('   ERROR: Invalid Path/File... Please make sure there are no spaces / special character(s) \n');
    fprintf('     in the filename. Press any key to exit...'); pause;
    fprintf('\n\n');
    return
end

try
    DRM(1) = DRM_DABNet(dt,A,B,NN,u_mean,y_mean,u_sigma,y_sigma,input_indices,output_indices,input_descs,output_descs,input_units,output_units,dt_unit);
catch
    fprintf('   ERROR: DABNet-Parameter file empty or incompatible type. If this is in error, please \n');
    fprintf('     see documentation for reporting-bugs. Press any key to exit...'); pause;
    fprintf('\n\n');
    return;
end
fprintf('\n');
fprintf('Conversion successful...written to workspace variable ''DRM'' \n');
fprintf('For preserving the DRM(s) for future use, please save the variable in a MAT file. \n\n');
clearvars -except DRM
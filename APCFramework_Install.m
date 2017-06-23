function APCFramework_Install
%% Installation File for CCSI APC Framework
%
% In order to run this tool, please run this file to setup the required
% paths. You MUST be in the current directory of this file!

% Code written and developed by:
% Priyadarshi Mahapatra, PhD
% URS Corporation / National Energy Technology Laboratory
% 04 June 2014

current_path = cd;
warning OFF BACKTRACE;

try
    cd('Toolbox');
catch
    error('You don''t appear to be in the APC Framework base directory');
end

fprintf('\n-------------------------------------------------------------------\n');
fprintf([' INSTALLING APC FRAMEWORK TOOLBOX R' ...
 sprintf('%1.2f',APCFrameworkVersion('release')) ...
 '(v' sprintf('%1.1f',APCFrameworkVersion('internal')) ')\n\n']);

if ~strcmpi(computer('arch'),'win32')
    warning('APCFramework:absent32bitArch',['32-bit MATLAB and/or WinOS not detected. The framework uses '...
             'Simulink - ACM/APD interface (distributed with Aspen Engineering Suite) '...
             'which is currently supported for 32-bit MATLAB only. You may still continue '...
             'with the installation but will be unable to use features / examples using this '...
             'interface for ACM/APD-based ''plant'' simulations.']);     
    fprintf('\n');
    in = input(' - Do you want to continue installation? (Not Recommended) (y/n): ','s');
    fprintf('\n');
    if ~(strcmpi(in,'y'))
        fprintf(' Installation Aborted. \n\n');
        cd(current_path);
        return;
    end
end
fprintf(' - Adding APC Framework Toolbox to MATLAB Search Path...');

generated_path = genpath(cd);
addpath(generated_path);

rehash
fprintf('Done\n\n');

in = input(' - Would you like to save the path changes? (Recommended) (y/n): ','s');
fprintf('\n');
if(strcmpi(in,'y'))
   warning off MATLAB:SavePath:PathNotSaved
   flagSaveFailed = savepath;
   if flagSaveFailed        
        warning('APCFramework:absentAdminRights',['It appears you do not have administrator rights on your computer to save the MATLAB path. '...
				 'This installation will be valid for current MATLAB session only. In order to run '...
                 'APC Framework you will need to install it each time you wish to use it. To fix this please '...
                 'contact your system administrator to obtain administrative privileges.']);
        fprintf('\n');
   end
end

fprintf(' APC Framework Toolbox Installation COMPLETE !\n');
fprintf('-------------------------------------------------------------------\n\n');

cd(current_path);
end

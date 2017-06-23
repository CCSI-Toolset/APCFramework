% CCSI APC Framework
% MATLAB's "ident" to APC Framework's DRM Converter Tool

% Code written and developed by:
% Priyadarshi Mahapatra
% URS Corporation / National Energy Technology Laboratory / U.S. Department of Energy
% Morgantown, WV 26505
% 20 Nov 2013

clear all; clc; close all;

fprintf('\n');
fprintf(' MATLAB''s Ident to APC Framework compatible DRM converter \n');
fprintf('---------------------------------------------------------- \n');
fprintf('This script converts System Identification Toolbox''s linear state-space \n');
fprintf('model(s) into APC Framework''s DRM_SS object(s) for use in MPC algorithm. \n');
fprintf('\n');
fprintf('MAT-file containing the Ident data MUST contain the following variables - \n');
fprintf('\n');
fprintf('1. IdentMod          - Identified Model (idss object). If more than one models are \n');
fprintf('                       present, this variable should be cell array of idss objects \n');
fprintf('   - OR - \n');
fprintf('   A, B, C, D, etc.  - Identified Model (state-space matrices). If more than one models \n');
fprintf('                       are present, each of these variables should be cell array \n');
fprintf('2. um_idx and ym_idx - DRM variable(s) index information, to identify DRM I/O \n');
fprintf('                       port(s) among plant port(s), as u_drm = u_plant(um_idx) \n');
fprintf('3. Uss and Yss       - Nominal variable values to compute physical values from \n');
fprintf('                       DRM''s deviation values \n');
fprintf('\n');
file_unspecified = 1;
while file_unspecified
    fprintf('Press any key to browse the MAT file containing Ident data...'); pause;
    [FileName,PathName] = uigetfile('*.mat','Select the file containing Ident data');
    fprintf('\n');
    if FileName == 0
        fprintf('  Error: You must provide a MAT file to continue or press Ctrl-C to exit.\n');
    else
        file_unspecified = 0;
    end
end
load([PathName FileName]); clear FileName PathName file_unspecified;
fprintf('\n');

% Determine number of DRM models
if isa(um_idx,'cell')
    nDRM = length(um_idx);
else
    nDRM = 1;
end

for jj = 1:nDRM
    try
        if isa(IdentMod,'cell')
            DRM(jj) = DRM_SS(IdentMod{jj});
        else
            DRM(jj) = DRM_SS(IdentMod);
        end
    catch me
        try
            if isa(A,'cell') && isa(B,'cell') && isa(C,'cell') && isa(D,'cell') && isa(x0,'cell') && isa(Ts,'cell') && isa(K,'cell')
                DRM(jj) = DRM_SS(A{jj},B{jj},C{jj},D{jj},x0{jj},Ts{jj},K{jj});
            else
                DRM(jj) = DRM_SS(A,B,C,D,x0,Ts,K);
            end
        catch me
            fprintf('ERROR: Empty or incompatible model type. Please see documentation for \n');
            fprintf('  state-space model specification. Press any key to exit...'); pause;
            fprintf('\n\n');
            clear nDRM jj
            return;
        end
    end
end

try
    for jj = 1:nDRM
        if isa(um_idx,'cell')
            DRM(jj).u_idx = um_idx{jj};
        else
            DRM(jj).u_idx = um_idx;
        end
        if isa(ym_idx,'cell')
            DRM(jj).y_idx = ym_idx{jj};
        else
            DRM(jj).y_idx = ym_idx;
        end
        if isa(Uss,'cell')
            DRM(jj).uss = Uss{jj};
        else
            DRM(jj).uss = Uss;
        end
        if isa(Yss,'cell')
            DRM(jj).yss = Yss{jj};
        else
            DRM(jj).yss = Yss;
        end
        if exist('um_name','var')
            if isa(um_name{1},'cell')
                DRM(jj).u_name = um_name{jj};
            else
                DRM(jj).u_name = um_name;
            end
        end
        if exist('ym_name','var')
            if isa(ym_name{1},'cell')
                DRM(jj).y_name = ym_name{jj};
            else
                DRM(jj).y_name = ym_name;
            end
        end
        if exist('um_unit','var')
            if isa(um_unit{1},'cell')
                DRM(jj).u_unit = um_unit{jj};
            else
                DRM(jj).u_unit = um_unit;
            end
        end
        if exist('ym_unit','var')
            if isa(ym_unit{1},'cell')
                DRM(jj).y_unit = ym_unit{jj};
            else
                DRM(jj).y_unit = ym_unit;
            end
        end
        if exist('tm_unit','var')
            if isa(tm_unit,'cell')
                DRM(jj).t_unit = tm_unit{jj};
            else
                DRM(jj).t_unit = tm_unit;
            end
        end
    end
catch me
    fprintf('ERROR: One of the required DRM variable missing or incompatible. Please see \n');
    fprintf('  documentation for correct ''identified-model'' specification. Press any key to exit...'); pause;
    fprintf('\n\n');
    clear nDRM jj
    return;
end

fprintf('Conversion successful...written to workspace variable ''DRM'' \n');
fprintf('For preserving the DRM(s) for future use, please save the variable in a MAT file. \n\n');
clearvars -except DRM
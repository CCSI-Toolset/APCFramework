% CCSI APC Framework
% MATLAB's "ident" to State-Space Converter Tool

% Code written and developed by:
% Priyadarshi Mahapatra
% URS Corporation / National Energy Technology Laboratory / U.S. Department of Energy
% Morgantown, WV 26505
% 26 May 2013

clc;

fprintf('\n');
fprintf(' MATLAB''s Ident to State-Space converter \n');
fprintf('------------------------------------------ \n');
fprintf('This script converts System Identification Toolbox''s ''ident'' data into \n');
fprintf('state-space matrices (A, B, C, etc.) which can then be accessed through workspace. \n');
fprintf('\n');
fprintf('The variable named ''IdentMod'' must be loaded into workspace. \n');
in = input('Is the correct ''IdentMod'' already present in the workspace? (y/n): ','s');
if ~(strcmpi(in,'y'))
    ident_unspecified = 1;
    while ident_unspecified
        fprintf('Press any key to browse the MAT file containing IdentMod(s)...'); pause;
        fprintf('\n');
        [FileName,PathName] = uigetfile('*.mat','Select the data file containing IdentMod(s)');
        if FileName == 0
            fprintf('  ERROR: You must provide a MAT file to continue or press Ctrl-C to exit. \n');
        else
            ident_unspecified = 0;
        end
    end
    load([PathName FileName]);
    fprintf(['-- Loaded IdentMod file ''' FileName ''' into workspace -- \n']);
    clear FileName PathName ident_unspecified
end
clear in

try
    if isa(IdentMod,'cell')
        nDRM = length(IdentMod);
    else
        nDRM = 1;
    end
catch me
    fprintf('   ERROR: ''IdentMod'' variable not present in workspace. Press any key to exit...'); pause;
    fprintf('\n\n');
    return;
end

try
    for jj = 1:nDRM
        if isa(IdentMod,'cell')
            assert(isa(IdentMod{jj},'idss'));
        else
            assert(isa(IdentMod,'idss'));
        end
    end
    if isa(IdentMod,'cell')
        fprintf('-- Reading %d %s model(s) from workspace -- \n',nDRM,class(IdentMod{1}));
    else
        fprintf('-- Reading %d %s model(s) from workspace -- \n',nDRM,class(IdentMod));
    end
catch me
    fprintf('   ERROR: ''IdentMod'' empty or incompatible type. Please see documentation for \n');
    fprintf('     ''IdentMod'' specification. Press any key to exit...'); pause;
    fprintf('\n\n');
    clear nDRM jj
    return;
end

% A = cell(1,nDRM); B = cell(1,nDRM);
% C = cell(1,nDRM); D = cell(1,nDRM);
% x0 = cell(1,nDRM); Ts = cell(1,nDRM);
% K = cell(1,nDRM);
% um_name_ = cell(1,nDRM); um_unit_ = cell(1,nDRM);
% ym_name_ = cell(1,nDRM); ym_unit_ = cell(1,nDRM);
% tm_unit_ = cell(1,nDRM);

if isa(IdentMod,'cell')
    for jj = 1:nDRM   
        A{jj} = IdentMod{jj}.A; B{jj} = IdentMod{jj}.B;
        C{jj} = IdentMod{jj}.C; D{jj} = IdentMod{jj}.D;
        x0{jj} = IdentMod{jj}.x0;
        Ts{jj} = IdentMod{jj}.Ts;
        K{jj} = IdentMod{jj}.K;
        if ~exist('um_name','var')
            um_name{jj} = strtrim(IdentMod{jj}.InputName);
        end
        if ~exist('ym_name','var')
            ym_name{jj} = strtrim(IdentMod{jj}.OutputName);
        end
        if ~exist('um_unit','var')
            um_unit{jj} = strtrim(IdentMod{jj}.InputUnit);
        end
        if ~exist('ym_unit','var')
            ym_unit{jj} = strtrim(IdentMod{jj}.OutputUnit);
        end
        if ~exist('tm_unit','var')
            tm_unit{jj} = strtrim(IdentMod{jj}.TimeUnit);
        end
    end
else
    A = IdentMod.A; B = IdentMod.B;
    C = IdentMod.C; D = IdentMod.D;
    x0 = IdentMod.x0;
    Ts = IdentMod.Ts;
    K = IdentMod.K;
    if ~exist('um_name','var')
        um_name = strtrim(IdentMod.InputName);
    end
    if ~exist('ym_name','var')
        ym_name = strtrim(IdentMod.OutputName);
    end
    if ~exist('um_unit','var')
        um_unit = strtrim(IdentMod.InputUnit);
    end
    if ~exist('ym_unit','var')
        ym_unit = strtrim(IdentMod.OutputUnit);
    end
    if ~exist('tm_unit','var')
        tm_unit = strtrim(IdentMod.TimeUnit);
    end
end

fprintf('Conversion successful...written to workspace variables ''A'', ''B'', ''C'', ''D'', ''Ts'', ''x0'' and ''K'' \n');
fprintf('For preserving the state-space parameters for future use, please save the variables in a MAT file. \n\n');
clear nDRM um_name_ ym_name_ um_unit_ ym_unit_ tm_unit_ jj

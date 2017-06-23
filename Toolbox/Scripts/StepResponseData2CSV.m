clear all; clc;

fprintf('\n');
fprintf(' Step Response Data to CSV Tool \n');
fprintf('---------------------------------- \n');
fprintf('This script converts Step Response Data using scripts like GetStepTestData.m into \n');
fprintf('into CSV which is compatible with user-specified CSV within DRM Builder Tool. \n');
fprintf('\n');
file_unspecified = 1;
while file_unspecified
    fprintf('Press any key to browse the MAT file containing Step Response Data...'); pause;
    fprintf('\n');
    [FileName,PathName] = uigetfile({'*.mat','MAT-files (*.mat)'},'Select the data file containing Step Responses');
     if FileName == 0
        fprintf('  ERROR: You must provide a MAT file to continue or press Ctrl-C to exit. \n');
    else
        file_unspecified = 0;
    end
end
load([PathName FileName]);
fprintf(['-- Loaded file ''' FileName ''' into workspace -- \n\n']);

file_unspecified = 1;
while file_unspecified
    fprintf('Press any key to specify desired path/filename for CSV file...'); pause;
    fprintf('\n');
    [FileName,PathName] = uiputfile({'*.csv','Comma Separated Value files (*.csv)'},'Save As...');
     if FileName == 0
        fprintf('  ERROR: You must provide a CSV file to continue or press Ctrl-C to exit. \n');
    else
        file_unspecified = 0;
    end
end

try
    csvFileID = fopen([PathName FileName],'w+');
catch me
    fprintf('ERROR: Unable to write to specified file...Check folder permissions...\n\n');
    return;
end
fprintf('\n');

clear file_unspecified me

if isa(U,'cell')
    nDRM = length(U);
else
    nDRM = 1;
end

if nDRM > 1
    fprintf('-- %d transient data-sets were found. \n\n',nDRM);    
    Mod_unspecified = 1;
    while Mod_unspecified
        Mod_idx = input('Enter the desired data-set to be used: ');
        if isempty(Mod_idx)
            Mod_idx = 1;
            Mod_unspecified = 0;
        elseif length(Mod_idx) > 1 || Mod_idx < 1 || Mod_idx > nDRM
            fprintf('  ERROR: You must provide one entry from within 1 to %d.\n', nDRM);
        else
            Mod_unspecified = 0;
        end
    end
else
    Mod_idx = 1;
end
clear nDRM Mod_unspecified
fprintf('-- Reading Data Set %d -- \n\n',Mod_idx);

%% Row 1
% Sampling Time, Ts
if isa(Ts,'cell')
    dt = Ts{Mod_idx};
else
    dt = Ts;
end
% Time Unit (optional)
if exist('tm_unit','var')
    if isa(tm_unit,'cell')
        t_unit = tm_unit{Mod_idx};
    else
        t_unit = tm_unit;
    end
else
    t_unit = '';
end
fprintf(csvFileID,[num2str(dt) ',']);
fprintf(csvFileID,'%s\n',t_unit);
clear dt t_unit

%% Row 2
% Number of I/O
if isa(U,'cell')
    nu = size(U{Mod_idx},1);
else
    nu = size(U,1);
end
if isa(Y,'cell')
    ny = size(Y{Mod_idx},1);
else
    ny = size(Y,1);
end
fprintf(csvFileID,[num2str(nu) ',' num2str(ny) '\n']);

%% Row 3
% Input Names
if exist('um_name','var')
    io_name = um_name{Mod_idx}{1};
    for i = 2:nu
        io_name = [io_name ',' um_name{Mod_idx}{i}];
    end
else
    io_name = 'u1';
    for i = 2:nu
        io_name = [io_name ',u' num2str(i)];
    end    
end
% Output Names
if exist('ym_name','var')
    for i = 1:ny
        io_name = [io_name ',' ym_name{Mod_idx}{i}];
    end
else
    for i = 1:ny
        io_name = [io_name ',y' num2str(i)];
    end    
end
fprintf(csvFileID,'%s\n',io_name);
%% Row 4
fprintf(csvFileID,'%s\n',io_name);
%% Row 5
fprintf(csvFileID,'%s\n',io_name);
clear io_name i

%% Row 6
% Input Units
if exist('um_unit','var')
    io_unit = um_unit{Mod_idx}{1};
    for i = 2:nu
        io_unit = [io_unit ',' um_unit{Mod_idx}{i}];
    end
else
    io_unit = [];
    for i = 2:nu
        io_unit = [io_unit ','];
    end    
end
% Output Units
if exist('ym_unit','var')
    for i = 1:ny
        io_unit = [io_unit ',' ym_unit{Mod_idx}{i}];
    end
else
    for i = 1:ny
        io_unit = [io_unit ','];
    end    
end
fprintf(csvFileID,'%s\n',io_unit);
clear io_unit

%% Row 7 onwards...
if isa(U,'cell')
    Um = U{Mod_idx};
else
    Um = U;
end
if isa(Y,'cell')
    Ym = Y{Mod_idx};
else
    Ym = Y;
end
N = min(size(Um,2),size(Ym,2));
fprintf('-- Found %d sampling point(s) -- \n\n',N);
k1_unspecified = 1;
while k1_unspecified
    k1 = input('Enter the 1st point: ');
    if isempty(k1)
        k1 = 1; k1_unspecified = 0;
    elseif length(k1) > 1 || k1 < 1 || k1 > N
        fprintf('  ERROR: You must provide one entry from within 1 to %d.\n',N);
    else
        k1_unspecified = 0;
    end
end
fprintf('-- First sample selected as %d -- \n\n',k1);
k2_unspecified = 1;
while k2_unspecified
    k2 = input('Enter the 2nd point: ');
    if isempty(k2)
        k2 = N; k2_unspecified = 0;
    elseif length(k2) > 1 || k2 < 1 || k2 > N
        fprintf('  ERROR: You must provide one entry from within 1 to %d.\n',N);
    else
        k2_unspecified = 0;
    end
end
fprintf('-- Last sample selected as %d -- \n\n',k2);

fprintf('Provide desired ''input-hold'' formatting...\n');
fprintf('  1. y[k+1] = f(y[k],u[k]) or u[k]: y[k] --> y[k+1] (default)\n');
fprintf('  2. y[k] = f(y[k-1],u[k]) or u[k]: y[k-1] --> y[k] \n');
in_hold_unspecified = 1;
while in_hold_unspecified
    in_hold = input('Enter your selection (1 or 2): ');
    if isempty(in_hold)
        in_hold = 1; in_hold_unspecified = 0;
    elseif length(in_hold) > 1 || (in_hold ~= 1 && in_hold ~=2)
        fprintf('  ERROR: You must select between one of these options (or press RETURN for default).\n');
    else
        in_hold_unspecified = 0;
    end
end

if in_hold == 1 
    for k = k1:k2
        io_val = num2str(Um(1,k));
        for i = 2:nu
            io_val = [io_val ',' num2str(Um(i,k))];
        end
        for i = 1:ny
            io_val = [io_val ',' num2str(Ym(i,k))];
        end
        fprintf(csvFileID,[io_val '\n']);
    end
else
    % Repeat the 1st input
    if k1 > 1    
        io_val = num2str(Um(1,k1-1));
    else
        io_val = num2str(Um(1,k1));
    end
    for i = 2:nu
        if k1 > 1    
            io_val = [io_val ',' num2str(Um(i,k1-1))];
        else
            io_val = [io_val ',' num2str(Um(i,k1))];
        end
    end
    for i = 1:ny
        io_val = [io_val ',' num2str(Ym(i,k1))];
    end
    fprintf(csvFileID,[io_val '\n']);    
    for k = k1+1:k2
        io_val = num2str(Um(1,k-1));
        for i = 2:nu
            io_val = [io_val ',' num2str(Um(i,k-1))];
        end
        for i = 1:ny
            io_val = [io_val ',' num2str(Ym(i,k))];
        end
        fprintf(csvFileID,[io_val '\n']);
    end
end

clear io_val Um Ym N k1 k2 k1_unspecified k2_unspecified in_hold in_hold_unspecified
fprintf(['-- Writing to file ''' FileName ''' -- \n\n']);
fclose(csvFileID);
clear csvFileID FileName PathName Mod_idx nu ny i k ans
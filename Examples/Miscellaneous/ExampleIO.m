% Example of very complex system of I/O variables
% -- Plant --
% Inputs  : [1 2 3 4 5 6 7]
% Outputs : [1 2 3 4 5 6 7]

clear all;

um_idx{1} = [3 7 1 4];
ym_idx{1} = [2 5 4 3 7];

um_idx{2} = [4 3 5 1 6 7];
ym_idx{2} = [4 7 1 3 2];

um_idx{3} = [1 7 5 3 4 2];
ym_idx{3} = [6 5 1 4 7 3];

ModL = length(um_idx);

for jj = 1:ModL
    fprintf('   -------- DRM %d -------- \n',jj);
    fprintf(['     Inputs  - [' num2str(um_idx{jj}) '] \n']);
    fprintf(['     Outputs - [' num2str(ym_idx{jj}) '] \n']);
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
        fprintf('     Error: You must provide indices corresponding to the available input(s).\n');
    else
        uc_unspecified = 0;
    end
end
yc_unspecified = 1;
while yc_unspecified
    yc_specified_idx = input('   Enter the desired control output(s) in array format (e.g. [1 2]): ');
    if ~all(ismember(yc_specified_idx,ym_available_idx))
        fprintf('     Error: You must provide indices corresponding to the available output(s).\n');
    else
        yc_unspecified = 0;
    end
end
for jj = 1:ModL
    [~,uc_idx{jj},~] = intersect(um_idx{jj},uc_specified_idx);
    [~,yc_idx{jj},~] = intersect(ym_idx{jj},yc_specified_idx);
    lc_specified_idx = um_idx{jj};
    lc_specified_idx(uc_idx{jj}) = [];
    [~,lc_idx{jj},~] = intersect(um_idx{jj},lc_specified_idx);
    % Just to verify all control variable among plant's variables are the same...
    u_control_idx{jj} = um_idx{jj}(uc_idx{jj});
    y_control_idx{jj} = ym_idx{jj}(yc_idx{jj});
end
clear um_available_idx ym_available_idx uc_specified_idx yc_specified_idx lc_specified_idx uc_unspecified yc_unspecified jj
clear ModL



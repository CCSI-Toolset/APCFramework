function [Z,A_idx] = findArrayIntersection(A)
% Finds intersection of multiple arrays
% Usage: [Z,A_idx] = findArrayIntersection(A)
%   where,
%       A is a cell-array containing the arrays
%       Z provides the intersected array
%       A_idx provides a cell-array providing the index information of Z in A

% Code written and developed by:
% Priyadarshi Mahapatra, PhD
% URS Corporation / National Energy Technology Laboratory, 2013
%
% Last Modified: 2 February 2014

    if iscell(A)
        N = length(A);
        A_idx = cell(size(A));
        [Z,~,A_idx{1}] = intersect(A{1},A{1});
        % First find the intersected array
        for i = 2:N
            Z = intersect(Z,A{i});
        end
        % Now find the indices for each array
        for i = 1:N
            [~,A_idx{i},~] = intersect(A{i},Z);
        end  
    else
        [Z,~,~] = intersect(A,A);
        [~,A_idx,~] = intersect(A,Z);
    end
end
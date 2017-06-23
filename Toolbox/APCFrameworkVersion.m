function toolver = APCFrameworkVersion(varargin)
% APC Framework Version Information
% Return APC Framework software version. This document also contains version update
% information.

% Code written and developed by:
% Priyadarshi Mahapatra, PhD
% URS Corporation / National Energy Technology Laboratory, 2013
%
% Version Revision List -
% [2014-06-02] - internal (1.5), release (2014.10), config (0.2), simulation (0.2)
% [2015-06-15] - internal (1.6), release (2015.06), config (0.3), simulation (0.3)
% [2015-10-31] - internal (1.7), release (2015.10), config (0.4), simulation (0.4)
% Last Modified: 21 Oct 2015

if nargin == 0
	toolver = 1.7;
else
	toolname = lower(varargin{1});
    switch toolname
        case 'internal' %APC Framework Internal Version
            toolver = 1.7;
        case 'release' % APC Framework Release Version
            toolver = 2015.10;
        case 'config' % APCConfigTool Version
            toolver = 0.4;
		case 'simulation' % APCSimulationTool Version
			toolver = 0.4;
        case 'subversion'
            toolver = 0;
        case 'date'
            toolver = 'November 30, 2015';
        otherwise
            error('Unknown APC Framework tool name specified.');
    end
end


%% Internal Version Changes
% v1.7
% - IPOPT NLP optimizer added
% - APCSimulationToolGUI added

% v1.6 (6/15/2015)
% - Major implementation changes
% - DABNet hessian function added

% v1.5 (6/15/2014)
% - Re written as an object orientated system
% - DABNet gradient & constraint jacobian functions added

% v0.1 (6/15/2012)
% - Initial APC Framework idea, implemented as a collection of m files


end


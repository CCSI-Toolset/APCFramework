classdef APCHistorian < handle

% Class-definition for APC related plots
%
% Code written and developed by:
% Priyadarshi Mahapatra, PhD
% URS Corporation / National Energy Technology Laboratory
%
% Last Modified: June 2015

properties
    T;                          % sample-times array
    R;                          % setpoint array
    Y;                          % plant-output array
    U;                          % plant-input array
    t_unit;                     % time unit
    y_name;                    % plant-output variable name
    y_unit;                    % plant-output variable unit
    u_name;                    % plant-input variable name
    u_unit;                    % plant-input variable unit
    Ym;                         % control-model output array [cell for MMPC]
    YfSD;                       % filtered control-model output SD array [cell for MMPC]
    W;                          % model-weights (MMPC only)
    TimerDB;                    % computational time database
    Integral_Error = 0;         % sum of residuals over entire simulation horizon
    
    % Flags related to user-desired plots
    fPlotControlVars = false;
    fPlotAllVars = false;
    fPlotCompTime = false;
    fPlotOptimIter = false;
    fPlotModelWeights = false;
    
    fPlotAll_DisplayModel = true;
    fPlotAll_DisplayFilterSD = true;
    
    % Handle for plots
    hDynamicPlot;
    hControllerVarsPlot;
    hAllVarsPlot;
    hModelWeightsPlot;
    hCompCostPlot;
    hOptimIterPlot;
    
    % Variables related to Runtime Dynamic Plot
    fShowDynamicPlot = false;
    fScrollPlot = false;
    FrameDuration;                  % scroll-width (in time variable)
    fFixedFrameInterval = false;    % flag for forcing fixed scroll-width
    fShowIntegralError = false;     % flag for displaying integral error
    fShowControlCalcTime = false;   % flag for displaying control calculation cost
end % properties

properties (SetAccess = private)
    APC;                    % main APC object    
end % properties    

properties (SetAccess = private, GetAccess = private)
%     APC;                    % main APC object    
    nyc; nuc;               % no. of control variables
    ModL;                   % no. of control-models
    Mod_idx;                % active model index among all control models
    fs = 12;                % default font-size
    fsb = 14;               % default font-size (bold)
    fw = 'bold';            % default font-weight [normal, bold, light, demi]
    lw = 1.8;               % default line-width (normal)
    lwl = 1.2               % default line-width (light)
    lwb = 2.4;              % default line-width (bold)
    ModelColor = {[0 0 1]; [0 0.5 0]; [1 0 0]; [0 0.75 0.75]; [0.75 0 0.75]; [0.75 0.75 0]; [0.25 0.25 0.25]};
    % Various Graphic Handles
    pl_R; pl_Y; pl_U;
    ax_Y; ax_U;
    txt_IE; txt_CT;
    
    IE_idx = 0;            % index through which mean integral error has been evaluated
    
end % properties
    
methods

    function obj = APCHistorian(APC,y_name,u_name,y_unit,u_unit,t_unit,T,R,Y,U,Ym)
        switch nargin
            case 1 % APCHistorian(APC)
                obj.APC = APC;
            case 6 % APCHistorian(APC,y_name,u_name,y_unit,u_unit,t_unit)
                obj.APC = APC;
                obj.y_name = y_name;
                obj.u_name = u_name;
                obj.y_unit = y_unit;
                obj.u_unit = u_unit;
                obj.t_unit = t_unit;
            case 10 % APCHistorian(APC,y_name,u_name,y_unit,u_unit,t_unit,T,R,Y,U)
                obj.APC = APC;
                obj.y_name = y_name;
                obj.u_name = u_name;
                obj.y_unit = y_unit;
                obj.u_unit = u_unit;
                obj.t_unit = t_unit;
                obj.T = T;
                obj.R = R;
                obj.Y = Y;
                obj.U = U;
                obj.evalControllerError();
            case 11 % APCHistorian(APC,y_name,u_name,y_unit,u_unit,t_unit,T,R,Y,U,Ym)
                obj.APC = APC;
                obj.y_name = y_name;
                obj.u_name = u_name;
                obj.y_unit = y_unit;
                obj.u_unit = u_unit;
                obj.t_unit = t_unit;
                obj.T = T;
                obj.R = R;
                obj.Y = Y;
                obj.U = U;
                obj.Ym = Ym;
                obj.evalControllerError();
            otherwise
                disp('Incorrect APCHistorian arguments provided. See documentation for correct usage.');
        end % switch
        
        % Some Nomenclature (for frequently occuring variables)
        obj.nyc = length(obj.APC.y_control_idx);
        obj.nuc = length(obj.APC.u_control_idx);
        if isa(obj.APC,'APC_MMPC')
            obj.ModL = obj.APC.ModL;
            obj.Mod_idx = obj.APC.Mod_idx;
        else
            obj.ModL = 1;
            obj.Mod_idx = 1;
        end
        if isempty(obj.Ym)
            obj.Ym = cell(1,obj.ModL); 
        end
        if isempty(obj.YfSD)
            obj.YfSD = cell(1,obj.ModL);
        end
    end % function

    function displayPlot(obj,varargin)
        for iarg = 1:nargin-1
            switch lower(varargin{iarg})
                case {'control variables','control vars','controlvars','control','cv'}
                    obj.plotControllerVariables;
                case {'all variables','all vars','allvars','all','av'}
                    obj.plotAllVariables;
                case {'computational time','computation time','computational cost','computation cost','compcost','comptime','ct'}
                    obj.plotComputationTime;
                case {'optimization iteration','optim iter','optimiter','oi'}
                    obj.plotOptimIteration
                case {'model weight','model weights','modweight','modwt','mw'}
                    obj.plotModelWeights;
                otherwise
                    disp('Invalid User-Plot Specified...');
            end % switch
        end % for
    end % function
    
    function setupPlot(obj,varargin)
        for iarg = 1:2:nargin-1
            switch lower(varargin{iarg})
                case {'showmodel','show model','showdrm','show drm'}
                    obj.fPlotAll_DisplayModel = varargin{iarg+1};
                case {'showfiltersd','showkfsd'}
                    obj.fPlotAll_DisplayFilterSD = varargin{iarg+1};
                otherwise
                    disp('Invalid Plot Configuration Specified...');
            end % switch
        end % for
    end % function    
                    
    function saveResponse(obj,varargin)
        fTriggerPlot = false;
        for iarg = 1:2:nargin-1 % scan for odd-numbered entries
            switch lower(varargin{iarg})
                case 'time'
                    t = varargin{iarg+1};
                    obj.T = [obj.T t];
                    fTriggerPlot = true; % update dynamic plot when time is updated
                case 'setpoint'
                    r = varargin{iarg+1};
                    obj.R = [obj.R r];
                case 'output'
                    y = varargin{iarg+1};
                    obj.Y = [obj.Y y];
                case 'input'
                    u = varargin{iarg+1};
                    obj.U = [obj.U u];
                case 'model'
                    ym = varargin{iarg+1};
                    if isa(ym,'cell') % typically for MMPC
                        for jj = 1:obj.ModL
                            obj.Ym{jj} = [obj.Ym{jj} ym{jj}];
                        end
                    else
                        obj.Ym{1} = [obj.Ym{1} ym];
                    end
                case 'filtersd'
                    yfSD = varargin{iarg+1};
                    if isa(yfSD,'cell') % typically for MMPC
                        for jj = 1:obj.ModL
                            obj.YfSD{jj} = [obj.YfSD{jj} yfSD{jj}];
                        end
                    else
                        obj.YfSD{1} = [obj.YfSD{1} yfSD];
                    end
                case 'weight'
                    w = varargin{iarg+1};
                    obj.W = [obj.W w];
                otherwise
                    disp('Incorrect option(s) provided. See documentation for correct usage.');
            end % switch
        end % for
        if fTriggerPlot && obj.fShowDynamicPlot && ishandle(obj.hDynamicPlot)
            N = length(obj.T);
            if obj.fScrollPlot
                fAboveScrollSize = false;
                Tbeg = obj.T(1); Tend = obj.T(end);
                if (Tend-Tbeg) > obj.FrameDuration
                    Tbeg = Tend - obj.FrameDuration;
                    fAboveScrollSize = true;
                else
                    Tend = Tbeg + obj.FrameDuration;
                end
            end
            for i = 1:obj.nyc
                i_ = obj.APC.y_control_idx(i);
                if ~isempty(obj.R)
                    set(obj.pl_R(i),'XData',obj.T,'YData',obj.R(i,1:N))
                end
                set(obj.pl_Y(i),'XData',obj.T,'YData',obj.Y(i_,1:N))
                if obj.fScrollPlot && (obj.fFixedFrameInterval || fAboveScrollSize)
                    set(obj.ax_Y(i),'XLim',[Tbeg Tend]);
                end
            end
            for j = 1:obj.nuc
                j_ = obj.APC.u_control_idx(j);
                set(obj.pl_U(j),'XData',obj.T,'YData',obj.U(j_,1:N))
                if obj.fScrollPlot && (obj.fFixedFrameInterval || fAboveScrollSize)
                    set(obj.ax_U(j),'XLim',[Tbeg Tend]);
                end                
            end
            
            if obj.fShowIntegralError
                set(obj.txt_IE,'String',['Mean Integral Error = ' num2str(obj.evalControllerError())]);
            end
            if obj.fShowControlCalcTime
                try % prevent error if TimerDB.Control_Calculation is non-existent
                    set(obj.txt_CT,'String',['Control Calculation Time = ' num2str(obj.TimerDB.Control_Calculation(end)) ' sec']);
                catch
                end
            end
            drawnow
        end % for
    end % function

    function showDynPlot(obj,varargin)
        for iarg = 1:2:nargin-1 % scan for odd-numbered entries
            switch lower(varargin{iarg})
                case {'frameduration','scrollwidth'}
                    obj.FrameDuration = varargin{iarg+1};
                    if obj.FrameDuration > 0
                        obj.fScrollPlot = true;
                    else
                        obj.fScrollPlot = false;
                    end
                case {'fixedframeinterval','forcefixedscrollwidth'}
                    if varargin{iarg+1} ~= true
                        obj.fFixedFrameInterval = false;
                    else
                        obj.fFixedFrameInterval = true;
                    end
                case {'showintegralerror'}
                    if varargin{iarg+1} ~= true
                        obj.fShowIntegralError = false;
                    else
                        obj.fShowIntegralError = true;
                    end                    
                case {'showcontrolcalctime'}
                    if varargin{iarg+1} ~= true
                        obj.fShowControlCalcTime = false;
                    else
                        obj.fShowControlCalcTime = true;
                    end                    
                otherwise
                    disp('Incorrect option(s) provided. See documentation for correct usage.');
            end % switch
        end

        %figure('units','normalized','outerposition',[0 0 1 1]);
        %obj.hDynamicPlot = figure('units','normalized','outerposition',[1/6 1/6 2/3 2/3]);
        pos_offset = 1/100; % offset from the top-left corner
        obj.hDynamicPlot = figure('units','normalized','outerposition',[pos_offset 1/3-2*pos_offset 2/3 2/3]);
        %figure;
        titleString = 'Dynamic Plot';
        if ~isempty(obj.APC.Tag)
            titleString = ['[' obj.APC.Tag '] ' titleString];
        end
        set(obj.hDynamicPlot,'Name',titleString,'NumberTitle','off');
        
        sp_no = 0;
        for i = 1:obj.nyc
            sp_no = sp_no + 1;
            i_ = obj.APC.y_control_idx(i);
            subplot(2,max(obj.nyc,obj.nuc),sp_no); hold on; grid on;
            obj.pl_R(i) = stairs(NaN,NaN,'r--','Linewidth',obj.lw,'DisplayName','Setpoint');
            obj.pl_Y(i) = plot(NaN,NaN,'k','Linewidth',obj.lwb,'DisplayName','Plant');
            obj.ax_Y(i) = gca;
            set(gca,'FontSize',obj.fs,'FontWeight',obj.fw);
            % Set y-Axis Labels
            if isempty(obj.y_name)
                ylabel(['y_' num2str(i_)],'FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                if ~isempty(obj.y_unit{i_})
                    ylabel([obj.y_name{i_} ' [' obj.y_unit{i_} ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
                else
                    ylabel(obj.y_name{i_},'FontSize',obj.fsb,'FontWeight',obj.fw);
                end
            end
            % Set x-Axis Labels
            if isempty(obj.t_unit)
                xlabel('Time','FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                xlabel(['Time [' obj.t_unit ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
            end
        end
        sp_no = max(obj.nyc,obj.nuc);
        for j = 1:obj.nuc
            sp_no = sp_no + 1;
            j_ = obj.APC.u_control_idx(j);
            subplot(2,max(obj.nyc,obj.nuc),sp_no); hold on; grid on;
            obj.pl_U(j) = stairs(NaN,NaN,'k','Linewidth',obj.lwb);
            obj.ax_U(j) = gca;
            set(gca,'FontSize',obj.fs,'FontWeight',obj.fw);
             % Set y-Axis Labels
            if isempty(obj.u_name)
                ylabel(['u_' num2str(j_)],'FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                if ~isempty(obj.u_unit{j_})
                    ylabel([obj.u_name{j_} ' [' obj.u_unit{j_} ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
                else
                    ylabel(obj.u_name{j_},'FontSize',obj.fsb,'FontWeight',obj.fw);
                end                    
            end
            % Set x-Axis Labels
            if isempty(obj.t_unit)
                xlabel('Time','FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                xlabel(['Time [' obj.t_unit ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
            end
        end
        
        axes('Visible','off','Unit','normalized','Position',[0 0 1 1]);
        if obj.fShowIntegralError
            obj.txt_IE = text(0.99,0.99,'','FontSize',obj.fs,'FontWeight',obj.fw,'Visible','on','Unit','normalized','HorizontalAlignment','right','VerticalAlignment','top');
        end
        if obj.fShowControlCalcTime
            obj.txt_CT = text(0.01,0.99,'','FontSize',obj.fs,'FontWeight',obj.fw,'Visible','on','Unit','normalized','HorizontalAlignment','left','VerticalAlignment','top','Rotation',0);
        end
        
        obj.fShowDynamicPlot = true;
    end % function
    
    function MIE = evalControllerError(obj,fReset)
        if nargin > 1 && fReset
            pause;
            obj.IE_idx = 0;
            obj.Integral_Error = 0;
        end
        IE = 0; N = length(obj.T);
        if ~isempty(obj.R)
        for k = obj.IE_idx+1:N-1
            Ts = obj.T(k+1) - obj.T(k);
            residual = obj.R(:,k) - obj.Y(obj.APC.y_control_idx,k);
            IE = IE + residual'*obj.APC.wy*residual * Ts/2; % Backward difference
            residual = obj.R(:,k) - obj.Y(obj.APC.y_control_idx,k+1);
            IE = IE + residual'*obj.APC.wy*residual * Ts/2; % Forward difference
        end
        end
        obj.Integral_Error = obj.Integral_Error + IE;
        obj.IE_idx = N-1;
        MIE = obj.Integral_Error/(obj.T(N)-obj.T(1));
    end

    function plotControllerVariables(obj)
        nyp = size(obj.Y,1);
        nup = size(obj.U,1);
        
        if length(obj.y_name) ~= nyp || length(obj.y_unit) ~= nyp
            disp('Incorrect output variable names or units specified. Please correct this for proper labels.');
            disp('Using default output notations (y1, y2, etc.)...');
        end
        if length(obj.u_name) ~= nup || length(obj.u_unit) ~= nup
            disp('Incorrect input variable names or units specified. Please correct this for proper labels.');
            disp('Using default input notations (u1, u2, etc.)...');
        end        
        if isempty(obj.t_unit)
            disp('Please provide the time units for proper labels');
        end
        
        obj.hControllerVarsPlot = figure('units','normalized','outerposition',[0 0 1 1]);
        titleString = 'Historian: Control Variables';
        if ~isempty(obj.APC.Tag)
            titleString = ['[' obj.APC.Tag '] ' titleString];
        end
        set(obj.hControllerVarsPlot,'Name',titleString,'NumberTitle','off');

        t = obj.T;
        N = length(t);
        sp_no = 0;        
        for i = 1:obj.nyc
            sp_no = sp_no + 1;
            i_ = obj.APC.y_control_idx(i);
            subplot(2,max(obj.nyc,obj.nuc),sp_no); hold all;
            if ~isempty(obj.R)
                r = obj.R(i,1:N);
                stairs(t,r,'r--','Linewidth',obj.lw,'DisplayName','Setpoint');
            end
            y = obj.Y(i_,1:N);
            plot(t,y,'k','Linewidth',obj.lwb,'DisplayName','Plant');
            if i == obj.nyc % -- Show legend only for last CV-subplot
                legend('-DynamicLegend','Location','Best');
            end
            % Set Axis Limits
            xlim([t(1) t(end)]); grid on;
            % Set global font (ticks, labels, title, etc.)
            set(gca,'FontSize',obj.fs,'FontWeight',obj.fw);
            % Set y-Axis Labels
            if length(obj.y_name) ~= nyp || length(obj.y_unit) ~= nyp
                ylabel(['y_' num2str(i_)],'FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                if ~isempty(obj.y_unit{i_})
                    ylabel([obj.y_name{i_} ' [' obj.y_unit{i_} ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
                else
                    ylabel(obj.y_name{i_},'FontSize',obj.fsb,'FontWeight',obj.fw);
                end
            end
            % Set x-Axis Labels
            if isempty(obj.t_unit)
                xlabel('Time','FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                xlabel(['Time [' obj.t_unit ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
            end
        end
        sp_no = max(obj.nyc,obj.nuc);
        for j = 1:obj.nuc
            sp_no = sp_no + 1;
            j_ = obj.APC.u_control_idx(j);
            u = obj.U(j_,1:N);
            subplot(2,max(obj.nyc,obj.nuc),sp_no);
            stairs(t,u,'k','Linewidth',obj.lwb);
            % Set Axis Limits
            y_range = max(u) - min(u);
            if y_range == 0
                y_range = abs(0.5*max(u));
            end
            ylim([min(u)-0.05*y_range max(u)+0.05*y_range]);
            xlim([t(1) t(end)]); grid on;
            % Set global font (ticks, labels, title, etc.)
            set(gca,'FontSize',obj.fs,'FontWeight',obj.fw);
            % Set y-Axis Labels
            if length(obj.u_name) ~= nup || length(obj.u_unit) ~= nup
                ylabel(['u_' num2str(j_)],'FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                if ~isempty(obj.u_unit{j_})
                    ylabel([obj.u_name{j_} ' [' obj.u_unit{j_} ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
                else
                    ylabel(obj.u_name{j_},'FontSize',obj.fsb,'FontWeight',obj.fw);
                end                    
            end
            % Set x-Axis Labels
            if isempty(obj.t_unit)
                xlabel('Time','FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                xlabel(['Time [' obj.t_unit ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
            end
        end
    end % function

    function plotAllVariables(obj,varargin)
        nyp = size(obj.Y,1);
        nup = size(obj.U,1);
        
        if length(obj.y_name) ~= nyp || length(obj.y_unit) ~= nyp
            disp('Incorrect output variable names or units specified. Please correct this for proper labels.');
            disp('Using default output notations (y1, y2, etc.)...');
        end
        if length(obj.u_name) ~= nup || length(obj.u_unit) ~= nup
            disp('Incorrect input variable names or units specified. Please correct this for proper labels.');
            disp('Using default input notations (u1, u2, etc.)...');
        end        
        if isempty(obj.t_unit)
            disp('Please provide the time units for proper labels');
        end
        
        u_idx = 1:nup;
        y_idx = 1:nyp;
        
        if nargin >= 2
            u_idx = varargin{1};
            if ~all(ismember(u_idx,1:nup))
                disp('Incorrect input plot indices provided. All input(s) will be plotted...');
                u_idx = 1:nup;
            end
        end
        if nargin >= 3
            y_idx = varargin{2};
            if ~all(ismember(y_idx,1:nyp))
                disp('Incorrect output plot indices provided. All output(s) will be plotted...');
                y_idx = 1:nyp;
            end                
        end
        d_idx = [];
        if nargin >= 4
            d_idx = varargin{3};
            if ~all(ismember(d_idx,1:nup))
                disp('Incorrect disturbance plot indices provided. Disturbance(s) will NOT be plotted separately...');
                d_idx = [];
            end
        end
            
        nu = length(u_idx); ny = length(y_idx); nd = length(d_idx);
        nRows = 2 + (nd > 0);
        
        % Find which plot to display legend
        i_legend = intersect(y_idx,obj.APC.y_control_idx);
        if isempty(i_legend)
            if ~isempty(y_idx)
                i_legend = y_idx(end);
            end
        else
            i_legend = i_legend(end);
        end
        
        obj.hAllVarsPlot = figure('units','normalized','outerposition',[0 0 1 1]);
        titleString = 'Historian: All Variables';
        if ~isempty(obj.APC.Tag)
            titleString = ['[' obj.APC.Tag '] ' titleString];
        end
        set(obj.hAllVarsPlot,'Name',titleString,'NumberTitle','off');
        
        t = obj.T;
        N = length(t);
        sp_no = 0;
        for i_plot = 1:ny
            i = y_idx(i_plot);
            sp_no = sp_no + 1;
            subplot(nRows,max([ny nu nd]),sp_no); hold all;
            if ~isempty(obj.R) && ismember(i,obj.APC.y_control_idx)
                i_ = i == obj.APC.y_control_idx;
                r = obj.R(i_,1:N);
                stairs(t,r,'r--','Linewidth',obj.lwl,'DisplayName','Setpoint');
            end
            for jj = 1:obj.ModL
                ii = obj.Mod_idx(jj);
                if obj.fPlotAll_DisplayModel && ~isempty(obj.Ym{jj}) && ismember(i,obj.APC.DRM(ii).y_idx)
                    i_ = i == obj.APC.DRM(ii).y_idx;
                    color_idx = rem(jj-1,length(obj.ModelColor))+1;
                    if isa(obj.APC,'APC_MMPC')
                        ModelName = ['Model-' num2str(ii)];
                    else
                        ModelName = 'D-RM';
                    end
                    ym = obj.Ym{jj}(i_,1:N);
                    if obj.fPlotAll_DisplayFilterSD && ~isempty(obj.YfSD{jj})
                        yfSD = obj.YfSD{jj}(i_,1:N);
                        hLower = plot(t,ym-yfSD,':','Linewidth',obj.lwl/2,'Color',obj.ModelColor{color_idx});
                        set(get(get(hLower,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
                        hUpper = plot(t,ym+yfSD,':','Linewidth',obj.lwl/2,'Color',obj.ModelColor{color_idx});
                        set(get(get(hUpper,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
                    end
                    plot(t,ym,'Linewidth',obj.lwl,'Color',obj.ModelColor{color_idx},'DisplayName',ModelName);
                end
            end
            y = obj.Y(i,1:N);
            plot(t,y,'k','Linewidth',obj.lw,'DisplayName','Plant');
%             if i == obj.APC.y_control_idx(end) % -- Show legend only for last CV-subplot
            if i == i_legend
                legend('-DynamicLegend','Location','Best');
            end
            % Set Axis Limits
            xlim([t(1) t(end)]); grid on;
            % Set global font (ticks, labels, title, etc.)
            set(gca,'FontSize',obj.fs,'FontWeight',obj.fw);
            % Set y-Axis Labels
            if length(obj.y_name) ~= nyp || length(obj.y_unit) ~= nyp
                ylabel(['y_' num2str(i)],'FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                if ~isempty(obj.y_unit{i})
                    ylabel([obj.y_name{i} ' [' obj.y_unit{i} ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
                else
                    ylabel(obj.y_name{i},'FontSize',obj.fsb,'FontWeight',obj.fw);
                end
            end
            % Set x-Axis Labels
            if isempty(obj.t_unit)
                xlabel('Time','FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                xlabel(['Time [' obj.t_unit ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
            end
        end
        
        sp_no = max([ny nu nd]);
        for j_plot = 1:nu
            j = u_idx(j_plot);
            sp_no = sp_no + 1;
            subplot(nRows,max([ny nu nd]),sp_no);
            u = obj.U(j,1:N);
            stairs(t,u,'k','Linewidth',obj.lw);
            % Set Axis Limits
            y_range = max(u) - min(u);
            if y_range == 0
                y_range = abs(0.5*max(u));
            end
            ylim([min(u)-0.05*y_range max(u)+0.05*y_range]);
            xlim([t(1) t(end)]); grid on;
            % Set global font (ticks, labels, title, etc.)
            set(gca,'FontSize',obj.fs,'FontWeight',obj.fw);            
            % Set y-Axis Labels
            if length(obj.u_name) ~= nup || length(obj.u_unit) ~= nup
                ylabel(['u_' num2str(j)],'FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                if ~isempty(obj.u_unit{j})
                    ylabel([obj.u_name{j} ' [' obj.u_unit{j} ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
                else
                    ylabel(obj.u_name{j},'FontSize',obj.fsb,'FontWeight',obj.fw);
                end                    
            end
            % Set x-Axis Labels
            if isempty(obj.t_unit)
                xlabel('Time','FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                xlabel(['Time [' obj.t_unit ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
            end
        end
        
      if nd > 0
        sp_no = 2*max([ny nu nd]);
        for j_plot = 1:nd
            j = d_idx(j_plot);
            sp_no = sp_no + 1;
            subplot(nRows,max([ny nu nd]),sp_no);
            u = obj.U(j,1:N);
            stairs(t,u,'k','Linewidth',obj.lw);
            % Set Axis Limits
            y_range = max(u) - min(u);
            if y_range == 0
                y_range = abs(0.5*max(u));
            end
            ylim([min(u)-0.05*y_range max(u)+0.05*y_range]);
            xlim([t(1) t(end)]); grid on;
            % Set global font (ticks, labels, title, etc.)
            set(gca,'FontSize',obj.fs,'FontWeight',obj.fw);            
            % Set y-Axis Labels
            if length(obj.u_name) ~= nup || length(obj.u_unit) ~= nup
                ylabel(['u_' num2str(j)],'FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                if ~isempty(obj.u_unit{j})
                    ylabel([obj.u_name{j} ' [' obj.u_unit{j} ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
                else
                    ylabel(obj.u_name{j},'FontSize',obj.fsb,'FontWeight',obj.fw);
                end                    
            end
            % Set x-Axis Labels
            if isempty(obj.t_unit)
                xlabel('Time','FontSize',obj.fsb,'FontWeight',obj.fw);
            else
                xlabel(['Time [' obj.t_unit ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
            end
        end
      end % if nd > 0
        
    end % function
    
    function plotModelWeights(obj,varargin)
        if nargin > 1
            obj.W = varargin{1};
        end
        try
            if isa(obj.APC,'APC_MMPC')
                t = obj.T;
                N = length(t);
                
                obj.hModelWeightsPlot = figure; hold all;
                titleString = 'Historian: Model Weights (MMPC)';
                if ~isempty(obj.APC.Tag)
                    titleString = ['[' obj.APC.Tag '] ' titleString];
                end
                set(obj.hModelWeightsPlot,'Name',titleString,'NumberTitle','off');
                
                for jj = 1:obj.ModL
                    ii = obj.Mod_idx(jj);
                    color_idx = rem(jj-1,length(obj.ModelColor))+1;
                    ModelName = ['Model-' num2str(ii)];
                    plot(t,obj.W(jj,1:N),'Linewidth',obj.lw,'Color',obj.ModelColor{color_idx},'DisplayName',ModelName);
                end
                legend('-DynamicLegend','Location','Best'); 
                ylim([0 1]);
                xlim([t(1) t(end)]); grid on;
                % Set global font (ticks, labels, title, etc.)
                set(gca,'FontSize',obj.fs,'FontWeight',obj.fw);
                % Set y-Axis Labels
                ylabel('Model Weight(s)','FontSize',obj.fsb,'FontWeight',obj.fw);
                % Set x-Axis Labels
                if isempty(obj.t_unit)
                    xlabel('Time','FontSize',obj.fsb,'FontWeight',obj.fw);
                else
                    xlabel(['Time [' obj.t_unit ']'],'FontSize',obj.fsb,'FontWeight',obj.fw);
                end
            end
        catch
            warning OFF BACKTRACE; warning('APCFramework:absentPlotModelWeights', ...
                ['Cannot plot model-weights. ''W'' non-existent / incompatible. '...
                'Please provide model-weight data (W) as an argument '...
                'to this plot function through <APCHistorian object>.plotModelWeights(W) or '...
                'supply data for current-step to <APCHistorian object>.W in runtime '...
                'through <APCHistorian object>.saveResponse(''Weight'',APC.W_model)']);
        end
    end % function
    
    function plotComputationTime(obj,varargin)
        if nargin > 1
            obj.TimerDB = varargin{1};
        end 
        try
            N = length(obj.TimerDB.Control_Calculation);

            obj.hCompCostPlot = figure; hold all;
            titleString = 'Historian: Computational Timers';
            if ~isempty(obj.APC.Tag)
                titleString = ['[' obj.APC.Tag '] ' titleString];
            end
            set(obj.hCompCostPlot,'Name',titleString,'NumberTitle','off');
                
            [AX,H1,H2] = plotyy(1:N,obj.TimerDB.Control_Calculation,1:N,obj.TimerDB.Plant_Simulation,'stairs','stairs');
            set(H1,'Linewidth',obj.lw); set(H2,'Linewidth',obj.lw);
            xlim(AX(1),[1 N]); grid on;
            xlim(AX(2),[1 N]);
            % Note: 'gca' will only affect LHS plot, using 'AX' instead
            set(AX,'FontSize',obj.fs,'FontWeight',obj.fw);
            title('Computation Cost','FontSize',obj.fsb,'FontWeight',obj.fw);
            xlabel('Time Step (k)','FontSize',obj.fsb,'FontWeight',obj.fw);
            ylabel(AX(1),'Control Calc. [sec]','FontSize',obj.fsb,'FontWeight',obj.fw);
            ylabel(AX(2),'Plant Sim. [sec]','FontSize',obj.fsb,'FontWeight',obj.fw);

%             figure; hold all;
%             stairs(1:N,obj.TimerDB.Control_Calculation,'Linewidth',obj.lw);
%             xlim([1 N]); grid on;
%             set(gca,'FontSize',obj.fs,'FontWeight',obj.fw);            
%             title('Computation Cost','FontSize',obj.fsb,'FontWeight',obj.fw);
%             xlabel('Time Step (k)','FontSize',obj.fsb,'FontWeight',obj.fw);
%             ylabel('Control Calc. [sec]','FontSize',obj.fsb,'FontWeight',obj.fw);

        catch
            warning OFF BACKTRACE; warning('APCFramework:absentPlotCompTimeData', ...
                ['Cannot plot computational-cost. ''TimerDB.Control-Calculation'' non-existent / incompatible. '...
                'Please provide computation-timer database (TimerDB) with valid Control-Calculation data as argument '...
                'to this plot function through <APCHistorian object>.plotComputationTime(TimerDB) or '...
                'supply data for current-step to <APCHistorian object>.TimerDB.Control-Calculation(k) in runtime.']);
        end
    end % function

    function plotOptimIteration(obj,varargin)
        if nargin > 1
            obj.TimerDB = varargin{1};
        end 
        try
            N = length(obj.TimerDB.Optim_Iteration);

            obj.hOptimIterPlot = figure; hold all;
            titleString = 'Historian: Optimization Performance';
            if ~isempty(obj.APC.Tag)
                titleString = ['[' obj.APC.Tag '] ' titleString];
            end
            set(obj.hOptimIterPlot,'Name',titleString,'NumberTitle','off');
                
            stairs(1:N,obj.TimerDB.Optim_Iteration,'Linewidth',obj.lw);
            xlim([1 N]); grid on;
            set(gca,'FontSize',obj.fs,'FontWeight',obj.fw);            
            title('Optimization Performance','FontSize',obj.fsb,'FontWeight',obj.fw);
            xlabel('Time Step (k)','FontSize',obj.fsb,'FontWeight',obj.fw);
            ylabel('# of iterations','FontSize',obj.fsb,'FontWeight',obj.fw);
        catch
            warning OFF BACKTRACE; warning('APCFramework:absentPlotOptimIterData', ...
                ['Cannot plot optimization-iteration. ''TimerDB.Optim_Iteration'' non-existent / incompatible. '...
                'Please provide computation-timer database (TimerDB) with valid Optim-Iteration data as argument '...
                'to this plot function through <APCHistorian object>.plotOptimIteration(TimerDB) or '...
                'supply data for current-step to <APCHistorian object>.TimerDB.Optim_Iteration(k) in runtime.']);
        end
    end % function
        
end % methods
    
end % classdef

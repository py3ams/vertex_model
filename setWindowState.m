function setWindowState(h,state)
% SETWINDOWSTATE sets a figures state to maximimze,minimize or restore
%
% SETWINDOWSTATE(H,STATE): H is figure handle or a vector of figure handles
%                          STATE is a string or cell array of strings-
%                               'maximize' - maximize figure
%                               'minimize' - minimize figure
%                               'restore'- restore figure
%                           if STATE is a string the state is applied to
%                           all H. If state is a cell array the length STATE
%                           must equal that of H and each state is applied
%                           individually.
%  Examples:
%   h= figure;
%   s = 'maximize';
%   setWindowState(h,s) %sets h to maximize
%
%   h(1) = figure;
%   h(2) = figure;
%   s = 'minimize';
%   setWindowState(h,s) %minimizes both figures
%
%   h(1) = figure;
%   h(2) = figure;
%   s = {'minimize','maximize'};
%   setWindowState(h,s) %minimizes h(1) and maximizes h(2)
%   setWindowState(h,'restore') %restores both windows
% Notes:
% 1. Figures must have 'Visible' set to 'on' and not be docked for
%    setWindowState to work.
% 2. Routine does not work for releases prior to R14SP2
% 3. The Java calls are undocumented by Mathworks
%
%Revisions: 01/09/06- Call the methods on the event thread using awtinvoke
%           05/11/06- Revisions for R2006a
%           09/28/06- Updated for R2006b and corrected warning call

warning off all

drawnow; %need to make sure that the figures have been rendered or Java error can occur

%check input argument number
error(nargchk(2, 2, nargin, 'struct'));

%is JVM available
if ~usejava('jvm')
	error('setWindowState requires Java to run.');
end

[j,s] = parseInput; %get the javaframes and desired operations
resizeWindow; %do the resizing operation

	function [j,s] = parseInput
		% is h all figure handles
		if ~all(ishandle(h)) || ~isequal(length(h),length(findobj(h,'flat','Type','figure')))
			error('All input handles must be valid figure handles');
		end %if
		
		%handle state argument
		if ischar(state)
			%make it a cell
			s = cellstr(repmat(state,[length(h),1]));
			
		elseif iscellstr(state)
			if length(state) ~= length(h)
				error('Cell array of strings: state must be same length as figure handle input');
			end %if
			s = state;
		else
			error('state must be a character array or a cell array of strings');
		end %if
		
		%check that the states are all valid
		if ~all(ismember(s,{'maximize','minimize','restore'}))
			error('Invalid states entered')
		end %if
		
		if length(h) == 1
			j{1} = get(h,'javaframe');
		else
			j = get(h,'javaframe');
		end %if
		
	end %parseInput

	function resizeWindow
		%get version so we know which method to call
		v = ver('matlab');
		%anticipating here that Mathworks will continue to change these
		%undocumented calls
		switch v(1).Release
			case {'(R14SP2)','(R14SP3)'}
				resize_method = 1;
			case {'(R2006a)','(R2006b)'}
				resize_method = 2;
			otherwise %warn but try method 2
				warning('setWindowState:UntestedRelease',['setWindowState has not been tested with release: ',v.Release]);
				resize_method = 2;
		end %switch
		
		switch resize_method
			case 1  %R14SP2-3
				for n = 1:length(j)
					awtinvoke(j{n},s{n});
				end %for
				
			case 2 %R2006a+
				for n = 1:length(j)
					fp= j{n}.fFigureClient.getFrameProxy;
					switch s{n}
						case 'maximize'
							awtinvoke(fp,'setMaximized(Z)',true)
						case 'minimize'
							awtinvoke(fp,'setMinimized(Z)',true)
						case 'restore'
							awtinvoke(fp,'setMaximized(Z)',false)
					end %switch
				end %for
			otherwise %should not happen
				error('Invalid resize method');
		end %switch
	end %resizeWindow

end %setWindowState
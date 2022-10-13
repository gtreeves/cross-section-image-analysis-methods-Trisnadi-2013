function soln = run_analyze_xs(filename,channels,channelnames,genotype,varargin)
%Runs "analyze_xs" on each lsm file in dir given in "pth".
%
%function soln = run_analyze_xs(filename,channels,channelnames,genotype,varargin)
%
% Optional argument varargin can consist of these things, in this order:
%	* "outfilename": If you want the output save files to have a name other
%		than just "run_analyze_xsCurrent.mat", it appends this 
%		character variable onto the end of the name.  Must be a character
%		variable, otherwise defaults.  Default, empty string, ''.  
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "fitprofiles": Logical variable for if you want the program to also
%		fit the secondary data (gene expression profiles) to canonical
%		profiles.   "true" means do the fitting, "false" means to skip the
%		fitting.  It is advantageous to skip the fitting if midlines need
%		to be found by hand and/or canonical peaks must be generated.
%		Default, "true".
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "yesplot": whether you want to plot the outcome.  Default, "false".  
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "ring_width": width of the annular ring.  That is, the distance in
%		microns into the embryo we take our data from. 
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "stage": what nuclear cycle is the embryo? it can be 10-14. Default,
%		14.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "nt": choice for number of bins in theta when detecting the embryo
%		periphery.  Default, 60. 
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "ns": choice for number of bins in "s" (pseudoarclength) when
%		measuring the intensity around the periphery.  Must be an even
%		number.  Default, 300.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.


%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	outfilename = varargin{iArg}; else
 	outfilename = '';
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	fitprofiles = varargin{iArg}; else
 	fitprofiles = true;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	yesplot = varargin{iArg}; else
 	yesplot = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	ring_width = varargin{iArg}; else
	ring_width = 18.36;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	stage = varargin{iArg}; else
 	stage = 14;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	nt = varargin{iArg}; else
	nt = 60;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	ns = varargin{iArg}; else
 	ns = 300;
end%, iArg = iArg + 1;

if ~ischar(outfilename)
	outfilename = '';
end
if ~(isnumeric(fitprofiles) || islogical(fitprofiles)) || ...
		~isscalar(fitprofiles) % must be numeric or logical, 1x1 vbl
	fitprofiles = true;
end

if ischar(filename) && ~isdir(filename)
	filenames = {filename};
elseif iscellstr(filename)
	filenames = filename;
elseif isdir(filename)
	pth = filename;
	filenames = readdir2(pth,'lsm');
else
	%
end
nLSM = length(filenames);

%
% Looping through and extract data from each file.
%
j = 1; LE = {}; FE = {}; e = 1;
for i = 1:nLSM
	filename = filenames{i};
	
	try
		data = analyze_xs(filename,channels,channelnames,genotype,...
			yesplot,ring_width,stage,nt,ns);
		if fitprofiles
			
			if any(channels == 3)
				data = fit_gaussian(data);
			end
			if any(channels == 2) || any(channels == 4)
				data = fit_peaks(data,false);
			end
		end
		save([filename(1:end-4),'_data.mat'],'data');
		soln(j,1) = data;
		save(['run_analyze_xsCurrent',outfilename],'soln')
		disp(['j = ',num2str(j)])
		j = j + 1;
	catch lastE
		LE{end+1} = lastE;
		FE{end+1} = filename;
		fprintf('%s in %s\n',LE{end}.message,filename)
		save(['run_analyze_xsError',outfilename],'LE','FE')
		disp(['e = ',num2str(e)])
		e = e + 1;
	end

end



if ~exist('soln','var')
	soln = 'no successful runs completed';
end







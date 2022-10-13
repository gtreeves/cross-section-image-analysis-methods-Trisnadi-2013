function [filenames,filenames_short] = readdir2(pth,ext)
%Outputs list of filenames found in "pth" with extension "ext".
%
%function [filenames,filenames_short] = readdir2(pth,ext)
%
% This function is "2" because I replaced using "strtok" with "deblank" and
% "strfind" on the "ispc" to handle the case where the filename has spaces
% in it.
%
% Inputs:
% 
% "pth": path to the directory you want (character variable)
% "ext": optional; character variable of extension of files you want.  If
%	you don't give it this variable, the default extension is "lsm".  
%	
%	NOTE: If you give '-' as the value of "ext", then the function will
%	look for files (or folders) without an extension (or, more precisely,
%	for files or folders that do not have any dots in their names). This
%	option currently only supported on windows machines.
%
% Outputs:
% 
% Both "filenames" and "filenames_short" are cell arrays that contain the
% files in your folder "pth" with extension "ext".  "filenames" includes
% the full path, and "filenames_short" only the filename itself.

if ~exist('ext','var')
	ext = 'lsm'; % default extension is "lsm"
end
if ~strcmp(ext(1),'.') && ~strcmp(ext,'-')
	ext = ['.',ext];
end

filenames = {};
filenames_short = {};
if ispc && ischar(pth)
	if ~strcmp(pth(end),char(92)) % if "pth" does not end in a backslash
		pth(end+1) = char(92); % then add the backslash.
	end
	w = ls(pth);

	%
	% Now we have a list of the filenames stored in "w".  This is an array
	% of size [m,n], where "m" is the number of entries in the directory,
	% and "n" is the number of characters in the longest filename (or
	% folder name).  The files that have shorter names have whitespace
	% appended afterward (to make "w" a rectangular array). We use a
	% function called strtok in a for loop that chops off the whitespace.
	% We will then store the separate strings in the cell vbl "filenames".
	%
	m = size(w,1);
	k = 1;
	for i = 1:m
		w1 = deblank(w(i,:));
		c1 = str2cell(w1,'.');
			
		if ~strcmp(ext,'-')
			name1 = cell2str(c1(1:end-1),'');
			ext1 = ['.',cell2str(c1(end),'')];
			record_criterion = strcmp(ext1,ext);
			
		else % no extension
			name1 = c1{1};
			ext1 = '';
			record_criterion = length(c1) == 1;
		end
		
		if record_criterion
			filenames{k,1} = [pth,name1,ext1];
			filenames_short{k,1} = [name1,ext1];
			k = k + 1;
		end
	end
	
elseif isunix && ischar(pth)
	if ~strcmp(pth(end),char(47)) % if "pth" does not end in a fwdslash
		pth(end+1) = char(47);
	end
	w = ls(pth);
	
	%
	% Now we have a list of the filenames stored in "w".  The files are
	% separated by whitespace, which can be either a tab (ASCII 9), space
	% (ASCII 32), or carriage return (ASCII 10 or 13?).  We use a function
	% called strtok in a while loop that parses our string for us.  We will
	% then store the separate strings in the cell vbl "filenames".
	%
	ww = w; k = 1;
	while ~isempty(ww)
		[tok,ww] = strtok(ww);
		% we have to filter this weird thing bc of strange unix behavior.
		kk = strfind(tok,char([27 91 48 48 109]));
		for kkk = length(kk):-1:1
			tok(kk(kkk):kk(kkk)+4) = [];
		end
		if length(tok) > length(ext)-1 && strcmp(tok(end-length(ext)+1:end),ext);
			filenames{k} = [pth,tok];
			filenames_short{k} = tok;
			k = k + 1;
		end
	end
end

filenames = sort(filenames);
filenames_short = sort(filenames_short);









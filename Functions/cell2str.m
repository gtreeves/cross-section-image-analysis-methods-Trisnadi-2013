function s = cell2str(c,sep)
%concatenates a row-cell array of strings into string vbl
%
%function s = cell2str(c,sep)
%
% This function takes a cell array, "c", and converts it into a character
% array, "s" that is essentially a comma (or other) separated list.  The
% elements of "c" get concatenated but separated by a special character,
% "sep".  For example, if "sep" is ',' (a comma), then "s" will become a
% comma separated list. If "c" is an m-by-n cell array, the substrings in
% the final output will be in the order of c(:) (long vector format).
%
% The sub-elements of "c" must all be single-row strings.  That is, if the
% character array of any element c{i} has multiple rows, this will not
% work as it is supposed to.
%
% This function is essentially the inverse of "str2cell".
%
% Ex: c = {'dog','cat','mouse'}
%	s = cell2str(c,',') returns s = 'dog,cat,mouse'.
%
% If "sep" is not passed, a comma is used as default.
%
% Copyright 2012, Dr Gregory T Reeves, Assistant Professor, NCSU

%
% Input checks
%
if ~iscellstr(c) && isempty(c)
	error('Input "c" must be a non-empty cell array of strings')
end
if ~exist('sep','var')
	sep = ',';
end
nsep = length(sep);

%
% Converting
%
c = c(:);
N = length(c);
for i = 1:N
	s1 = c{i};
	[m1,n1] = size(s1);
	if m1 > 1
		error('elements of input "c" must be single-row strings')
	end
	c{i} = [s1,sep]';
end
s = char(c)';

if nsep > 0
	s(end-nsep+1:end) = [];
end



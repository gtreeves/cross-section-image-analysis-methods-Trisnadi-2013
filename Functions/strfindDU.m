function v = strfindDU(cellstr,pattern)
%Just like "strfind" but gives true-false (rather than indices)
%
% This function also can take an N-by-1 cell array for "cellstr" and a
% 1-by-p cell array for "pattern" and spit out an N-by-p logical array,
% where the first column is like strfindDU(cellstr,pattern{1}), and the
% second column is like strfindDU(cellstr,pattern{2}), etc.

[m,n] = size(cellstr);
[r,p] = size(pattern);
if iscellstr(pattern) && r == 1 && n == 1
	v = false(m,p);
	for j = 1:p
		C = strfind(cellstr,pattern{j});
		for i = 1:m
			if ~isempty(C{i})
				v(i,j) = true;
			end
		end
	end
	
elseif ischar(pattern) && r == 1

	C = strfind(cellstr,pattern);
	v = false(m,n);
	for i = 1:m
		for j = 1:n
			if ~isempty(C{i,j})
				v(i,j) = true;
			end
		end
	end
elseif isempty(cellstr)
	v = [];
else
	error('You''ve screwed something up somewhere.')
end
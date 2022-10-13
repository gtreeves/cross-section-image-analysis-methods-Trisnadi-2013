function [x,A,B] = mDU(X,s)
%normalizes your n-by-1 or 1-by-n vector to fall b/w zero and one.
%
% function [x,A,B] = mDU(X,s)
%
% "s" could be a character variable, either 'row' or 'col' that specifies
% that the rows or columns of X should be done separately.  If nonexistent
% or not 'row' or 'col', then X will be done as a whole.
%
% x = (X - B)/(A - B); 
%	where A = max(X(:)) and B = min(X(:)), etc

if exist('s','var')
	N = size(X);
	x = zeros(N);
	if strmatch(s,'row')
		A = zeros(N(1),1); B = A;
		for i = 1:N(1)
			A(i) = max(X(i,:)); B(i) = min(X(i,:));
			x(i,:) = (X(i,:) - B(i))/(A(i) - B(i));
		end
	elseif strmatch(s,'col')
		A = zeros(1,N(2)); B = A;
		for i = 1:N(2)
			A(i) = max(X(:,i)); B(i) = min(X(:,i));
			x(:,i) = (X(:,i) - B(i))/(A(i) - B(i));
		end
	end
else
	A = max(X(:)); B = min(X(:));
	x = (X - B)/(A - B);
end
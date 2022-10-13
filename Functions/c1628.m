function varargout = c1628(I16,fh,Imax,Imin)
% takes an RGB uint16 and transforms it into uint8.
%
% function I8 = c1628(I16,fh,Imax,Imin)
%
% if no output is asked for, none is given, and you get an imshow instead.
%
% This function also works if "I16" is a double (ie, not just uint16).
%
% "fh" is a figure handle, if you want to send your plot to a particular
%	figure.  If this isn't specified, and no output is asked for, then
%	this function will just say "figure" rather than "figure(fh)".
%
% "Imax,Imin" are DOUBLE variables that, if specified, serve as the max and
% min for the whole image stack.  This is only important if you have an
% image stack of 2D images that all represent the same channel (but a
% z-stack or a time course).  Otherwise, it makes no sense to specify that.

I = double(I16);
clear I16
[m,n,o] = size(I);
if ~exist('Imin','var')
	nomaxmin = true;
else
	nomaxmin = false;
end
for i = 1:o
	I1 = I(:,:,i);
	if nomaxmin
		Imax = max(I1(:));
		Imin = min(I1(:));
	end
	
	I(:,:,i) = 255*(I1 - Imin)/(Imax - Imin);
end
I8 = uint8(I);


if nargout == 0
	if exist('fh','var') && ~isempty(fh) && round(fh) == fh
		if(any(findall(0,'Type','Figure') == fh)) % check if figure is already open
			s = get(fh,'visib');
			if strcmp(s,'off')
				figure(fh)
				set(fh,'visib','off')
			end
		else
			figure(fh)
		end
	else
		figure
	end
	if size(I8,3) == 2 % only two channels: red and green.
		Z = zeros(size(I8,1),size(I8,2));
		I8(:,:,3) = Z;
	elseif size(I8,3) > 3 % not grayscale or RGB, so sum all together.
		I8 = sum(double(I8),3);
		Imax = max(I8(:));
		Imin = min(I8(:));
		I8 = 255*(I8 - Imin)/(Imax - Imin);
	end
	imshow(I8,[])
else
	varargout = {I8};
end
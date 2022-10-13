function t_out = subtrbkgrnd(I,rbr,smoothened_already)
%Smooths I data, then subtracts background.
%
%function t_out = subtrbkgrnd(I,rbr,smoothened_already)
%
% This function re-processes the I data in case a different rbr is needed
% after the background had been subtracted.
%
% "I": can be a 1D image, a cell array of 1D images, or a multi-dimensional
%	structure, one field of which must be "Smth" (which is a cell array of
%	1D images).
%
% "t_out": 

if ~exist('smoothened_already','var')
	smoothened_already = true;
end

if isstruct(I)
	t_out = I;
	o = length(I);
	
	for k = 1:o
		if isfield(I,'metadata')
			t = I(k).t;
			t_out(k).t = sb(t,rbr,smoothened_already);
		else
			m = length(I(k).Smth); % "Smth" is a cell array.
			Smth = I(k).Smth;
			
			for j = 1:m
				n = length(Smth);
				I2 = cell(n,1);
				
				for i = 1:n
					I2{i} = sb(Smth{i},rbr,smoothened_already);
				end
			end
			t_out(k).t = I2;
			t_out(k).rbr = rbr;
		end
	end
	
elseif iscell(I)
	n = length(I);
	t_out = cell(n,1);
	
	for i = 1:n
		t_out{i} = sb(I{i},rbr,smoothened_already);
	end
	
elseif isnumeric(I)
	t_out = sb(I,rbr,smoothened_already);
else
	error('You screwed up.')
end



% -------------- subfunction to actually do the subtraction -----------
function t_out = sb(I,rbr,smoothened_already)

[ns,o] = size(I); ns = ns - 1;
rbr = rbr*ns; % "rolling ball radius"
t_out = zeros(ns+1,o);

for j = 1:o
	t1 = [I(1:end-1,j);I(:,j);I(2:end,j)];
	
	if ~smoothened_already
		t1 = smooth(I2,10);
	end
	
	se = strel('line',rbr(j),90);
	Itop = imtophat(t1,se);
% 	Iop = imopen(t1,se);
% 	Ie = imerode(t1,se);
% 	Itop = imsubtract(t1,Ie);
	t_out(:,j) = Itop(ns:2*ns);
end






function t = smooth_intron(S,T,s)
%Smooths out intronic probe data.
%
%function t = smooth_intron(S,T,s)
%
% This function smooths out intronic probe data.  The reason why is
% two-fold.  First, intronic data is inherently salt-and-pepper, so we are
% trying to convert this into something smooth.  Second, since intronic
% probe data (i.e., "nuclear dots") are one data point-per-nucleus, we are
% trying to convert these data into a smooth curve with points
% equally-spaced in "s".
%
% Inputs:
% "S": locations of all of the nuclei
% "T": intensity of intronic dots for each nucleus. If more than one
%	intrnoic probe is present (i.e., more than one color channel), then
%	there will be multiple columns of "T".
% "s": equally-spaced DV coordinate mesh from -1 to 1.
%
% Output:
% "t": the smoothed intronic data.

ns = length(s) - 1;
t1 = zeros(ns,1);
for i = 1:ns
	v = S > s(i) & S < s(i+1);
	if any(v)
		T1 = sort(T(v));
		n = min(length(T1),5);
		t1(i) = mean(T1(n:end));
	else
		t1(i) = NaN;
	end
end
p = 50;
t1 = [t1(end-p+1:end);t1;t1(1:p)];

%
% A very complicated way to replace NaN's with reasonable numbers.
%
k = find(isnan(t1));
for i = k(:)'
	k1 = max(1,i-1);
	while isnan(t1(k1))
		k1 = max(1,k1-1);
		if k1 == 1
			k1 = [];
			break
		end
	end
	k2 = min(length(t1),i+1);
	while isnan(t1(k2))
		k2 = min(length(t1),k2+1);
		if k2 == length(t1)
			k2 = [];
			break
		end
	end
	t1(i) = sum(t1([k1 k2]))/length([k1 k2]);
	
end

%
% Periodically extending "s" so we can smooth our data
%
s1 = s(1:end-1);
s1 = [s1(end-p+1:end)-2;s1;s1(1:p)+2];
tsmooth = smooth(s1,t1,floor(p/2));
t = tsmooth(p+1:end-p+1);






function [sV,sD,s_offset] = canonicalgeneborders(genename,genotype,hh)


matname = [genename,'avg'];
			
load(matname,'s','t','s_offset')
i0 = find(t == 0);
[imaxt,imaxt] = max(t);
i0V = i0(i0 < imaxt);
if ~isempty(i0V)
	i0V = i0V(end);
	sV = spline(t(i0V:imaxt),s(i0V:imaxt),hh);
else
	sV = NaN;
end
i0D = i0(i0 > imaxt);
if ~isempty(i0D)
	i0D = i0D(1);
	sD = spline(t(imaxt:i0D),s(imaxt:i0D),hh);
else
	sD = NaN;
end
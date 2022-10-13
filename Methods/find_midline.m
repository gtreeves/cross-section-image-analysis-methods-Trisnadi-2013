function s_mid = find_midline(s,T,S0,s_mid_in)
%Finds the ventral midline from gene expression data.
%
%function s_mid = find_midline(s,T,S0,s_mid_in)
%
% This function takes gene expression patterns and determines where the
% ventral midline should be.  The genes are grouped into three categories:
% dorsal genes, ventral genes, and lateral genes.  This function uses these
% categories to help determine the ventral midline.  Problems arise when
% there are more than one gene in a given color channel.  In that case, we
% turn to the channels that have single genes to help us out.
%
% "s": DV coordinate, from -1 to 1. "T": gene expression intensities in
% each relevant color channel. "S0": presumptive locations of each peak of
% gene expression.  If an
%	element of S0 = 1 or 0, then that means there is a gene in that color
%	channel that is eithe fully dorsal (such as dpp) or fully ventral
%	(e.g., sna).
% "s_mid_in": the present estimate of the midline.  If this is NaN, then
%	the program will continue with its estimate of the midline. Otherwise,
%	the program will skip this entry so as not to perturb the
%	previously-found midline.
%
% Output: "s_mid": the presumed location of the ventral midline for each
% color
%	channel. This is a 1-by-n vector where n is the number of color
%	channels that have mRNA in them.

n_ch = length(S0);
S_mid = cell(n_ch,1);
n_smid = zeros(1,n_ch);

%
% Looping through each channel that contains mRNA
%
for j = 1:n_ch
	if isnan(s_mid_in(j))
		s0 = S0{j};
		t = T(:,j);
		
		%
		% The number of peaks to look for, and the way we look for them,
		% depends on the types of genes.  Ventral or dorsal genes each get
		% one peak, while lateral get two. If we have the case of only
		% ventral or dorsal peaks, then our midline will be estimated by
		% the symmetry score.  If we have only one single lateral gene,
		% then we'll
		%
		n_genes = length(s0);
		v = s0 == 1 | s0 == 0;
		n_peaks = sum(v) + 2*sum(~v);
		
		%
		% Finding peak boundaries
		%
		h = 0.15; % noisethresh = h*max(t);
		sB = zeros(n_peaks,2);
		tmax = zeros(n_peaks,1);
		imax = tmax;
		for k = 1:n_peaks
			[tmax(k),imax(k)] = max(t);
			idx = find((t-h*tmax(k)).*(t([2:end,1])-h*tmax(k)) < 0);
			
			%
			% averaging to find peak bounds. (ie, not doing linear
			% interpolation.)  The boundaries are calculated at a cutoff of
			% peak height times "h".  We have to find these boundaries so
			% that we can set the value of the peak equal to zero in order
			% to find all the peaks.
			%
			if imax(k) < idx(1) || imax(k) > idx(end) % the peak is split
				
				sB(k,1) = 0.5*(s(idx(end)) + s(idx(end)+1));
				sB(k,2) = 0.5*(s(idx(1)) + s(idx(1)+1));
				
				t(idx(end):end) = 0;
				t(1:idx(1)+1) = 0;
			else
				jf = find(idx < imax(k));
				if isempty(jf)
					jf = 1;
				end
				jf = jf(end);
				
				sB(k,1) = 0.5*(s(idx(jf)) + s(idx(jf)+1));
				sB(k,2) = 0.5*(s(idx(jf+1)) + s(idx(jf+1)+1));
				
				if idx(jf+1) == length(t)
					t(idx(jf):end) = 0;
				else
					t(idx(jf):idx(jf+1)+1) = 0;
				end
			end
		end
		
		%
		% Now that we have the boundaries of all the peaks, we will average
		% them to find the peak locations.
		%
		sB(:,2) = sB(:,2) + 2*(sB(:,1) > sB(:,2));
		sP = sB*[0.5;0.5];
		sP = mod(sP+1,2) - 1; % making sure we don't split peaks.
		
		%
		% There are several different cases depending on what kind of genes
		% are in your channel(s).  Depending on what case we're in will
		% determine how we estimate s_mid from our data.
		%
		if all(s0 == 0)
			%
			% This is the simple case in which the ventral midline is just
			% the location of our peak of the ventral gene.  There should
			% only be one gene and one peak in this case.
			%
			S_mid{j} = sP;
			n_smid(j) = length(S_mid{j});
			
		elseif all(s0 == 1)
			%
			% This is the simple case in which the ventral midline is just
			% the complement of the location of our peak of the dorsal
			% gene.  There should only be one gene and one peak in this
			% case.
			%
			sDM = sP;
			S_mid{j} = mod(sDM,2)-1;
			n_smid(j) = length(S_mid{j});
			
		elseif sum(v) == n_genes
			%
			% This is the case in which there are two genes, one completely
			% ventral and the other dorsal. We don't have a way to be sure
			% which is which at this point, so we'll be ambivalent about it
			% at the moment.
			%
			S_mid{j} = sP;
			n_smid(j) = length(S_mid{j});
			
		elseif ~any(s0 == 0) && ~any(s0 == 1)
			%
			% This is the case in which we have just lateral genes.  In
			% this case, we assume there is only one gene in the channel,
			% because having two lateral genes in the same channel is
			% folly.  And because lateral genes are always closer to the
			% ventral midline than the dorsal midline, we take the ventral
			% midline as the average between the two peak locations.  We
			% have to take into account the provision that maybe the
			% ventral midline is on the other side of the graph from one of
			% the peaks.
			%
			if length(sP) ~= 2
				error('Gene channel with lateral gene only should have only 2 peaks')
			end
			dsP = abs(diff(sP));
			if dsP > 1 % then our peak is split.
				dsP = 2 - dsP;
				s_mid = max(sP) + 0.5*dsP;
			else % our peak is normal
				s_mid = min(sP) + 0.5*dsP;
			end
			S_mid{j} = mod(s_mid+1,2)-1;
			n_smid(j) = length(S_mid{j});
			
		elseif n_genes > 1 && any(s0 == 0)
			%
			% This is the case in which there are two genes, one completely
			% ventral and the other lateral.  In this case, the ventral
			% gene is the "middle" peak, and the lateral gene is the other
			% two peaks. We'll take the average between where the two genes
			% predict the vm to be.
			%
			% The ventral gene predicts the ventral midline to be at its
			% peak.
			%
			% The lateral gene /usually/ predicts the ventral midline to be
			% mid-way in between its two lateral peaks.  I say usually
			% because most of the time, the two lateral peaks are closer to
			% the vm than the dm.  There are some cases in which the two
			% peaks of the lateral gene are past 50% V-to-D.  In that case,
			% the mid-way calculation will give us instead the dorsal
			% midline rather than the ventral midline.  So we inquire
			% whether this is giving us the dm or the vm by which is closer
			% to the prediction from the ventral gene.
			%
			% We distinguish the two genes by how far three peaks are from
			% each other.  The ventral gene has a peak closer to the other
			% two peaks than the other two peaks are to each other.
			%
			
			%
			% Here we figure out how far the three peaks are all from each
			% other.  This will help us figure out which one is the ventral
			% gene.
			%
			sP1 = sort(sP);
			dsP = abs(diff(sP1([1:end,1])));
			v = dsP > 1;
			dsP(v) = 2 - dsP(v);
			
			%
			% The three cases.  If the biggest difference between the two
			% peaks is dsP(1) (so imax = 1), then the ventral gene is peak
			% 3.  If imax = 2, then the ventral gene is peak 1.  Else, the
			% ventral gene is peak 2.  This will tell us how the peaks are
			% split.
			%
			[dmax,imax] = max(dsP);
			switch imax
				case 1
					s_mid1 = sP1(3);
					if v(1)
						s_mid2 = mean([sP1(1)+2,sP1(2)]);
					else
						s_mid2 = mean(sP1(1:2));
					end
				case 2
					s_mid1 = sP1(1);
					if v(2)
						s_mid2 = mean([sP1(2),sP1(3)-2]);
					else
						s_mid2 = mean(sP1(2:3));
					end
				case 3
					s_mid1 = sP1(2);
					s_mid2 = mean(sP1([1 3]));
			end
			
			%
			% Now we ask whether the midling predicted by the lateral gene
			% is the dm or vm.  We'll figure this out by comparing it to
			% what the ventral gene predicts.
			%
			s_mid22 = mod(s_mid2,2)-1;
			if abs(s_mid1 - s_mid2) < abs(s_mid1 - s_mid22)
				S_mid{j} = mean([s_mid1,s_mid2]);
				n_smid(j) = length(S_mid{j});
			else
				S_mid{j} = mean([s_mid1,s_mid22]);
				n_smid(j) = length(S_mid{j});
			end
			
		elseif n_genes > 1 && any(s0 == 1)
			%
			% As far as I can tell, this is the final case. Here we have
			% two genes, one completely dorsal and the other lateral.  This
			% case is more difficult to figure out.  We will call each peak
			% a potential sDM (dorsal midline).  We will be ambivalent
			% about it after that.
			%
			sDM = sP;
			S_mid{j} = mod(sDM,2)-1;
			n_smid(j) = length(S_mid{j});
		end
	else
		S_mid{j} = s_mid_in(j);
		n_smid(j) = 1;
	end
end

s_mid = zeros(1,n_ch);
for j = find(n_smid == 1)
	s_mid(j) = S_mid{j};
end
s_star = mean(s_mid(n_smid == 1));
for j = find(n_smid ~= 1)
	svm = S_mid{j};
	[m,imin] = min(abs(svm - s_star));
	s_mid(j) = svm(imin);
end



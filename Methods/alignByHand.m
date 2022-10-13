function soln = alignByHand(soln,idx,rbr,plotintronic)
%aligns the midline of soln(idx) by hand.

if ~exist('idx','var') || isempty(idx)
	idx = 1:length(soln);
end
if ~exist('plotintronic','var')
	plotintronic = false;
end

for j = 1:length(idx)
	data = soln(idx(j));
	channels = data.channels;
	channels(channels == 5) = [];
	
	nmid = length(data.s_mid);
	if plotintronic
		%
		% Extracting nuclear (intronic) data
		%
		S = data.S;
		Intron = data.Intron;
		s = []; intron = [];
		for i = 1:length(Intron);
			s = [s;S{i}];
			intron = [intron;Intron{i}];
		end
		S = s;
	end
	
	
	%
	% Extracting smear data
	%
	T = data.t;
	T(:,channels == 1) = [];
	if exist('rbr','var') && length(rbr) == size(T,2)
		T = subtrbkgrnd(T,rbr);
	end
	s = data.s;	
	
	%
	% Dividing "s" in half
	%
	ns = length(s) - 1;
	if isodd(ns)
		error('you have an even number of points in x...eh?!?')
	end
	npts = ns/2 + 1;
	s1 = linspace(0,1,npts)';
	
	%
	% Initializing
	%
	figure('Position',[50 70 1500 500])
	fignum = gcf;
	button = 1;
	s_mid = 0;	
	while button == 1
		
		%
		% The normal plot
		%
		subplot(1,2,1)
		if plotintronic
			plot(S,intron,'. m')
			hold on			
		end
		h = plot(s,T,'Linewidth',2);
% 		set(h(1),'Color','b')
% 		set(h(2),'Color','r')
		hold off
		xlim([-1 1])
		title('Click over here')
		
		%
		% The split plot
		%
		t1 = T(npts:end,:); 
		t2 = flipud(T(1:npts,:));
		subplot(1,2,2)
		if plotintronic
			S1 = abs(S);
			h = plot(S1,intron,'.');
			set(h(1),'Color','c')
			hold on
		end
		h = plot(s1,[t1,t2],'Linewidth',2);
		set(h(1),'Color','b')
		set(h(2),'Color','r')
		set(h(3),'Color',[0 0.5 0])
		set(h(4),'Color','m')
		
		hold off
		xlim([0 1])
		title('Don''t click over here')
		
		%
		% User interaction
		%
		[s_star,y,button] = ginput(1);
		if button ~= 1
			break
		end
		
		%
		% Realigning each slice.
		%
		s_mid = s_mid + s_star;
		s_mid = mod(s_mid+1,2) - 1;
		s_comp = mod(s_star,2) - 1; % the complement.
		if plotintronic
			S = mod(S-s_star+1,2) - 1;
		end
		
		[i_comp,i_comp] = roundx(s_comp,s);
		T(end,:) = [];
		if i_comp == ns+1, i_comp = 1; end
		T = T([i_comp:end,1:i_comp],:);		
	end
	close(fignum)
	soln(idx(j)).s_mid = repmat(s_mid,1,nmid);	
end








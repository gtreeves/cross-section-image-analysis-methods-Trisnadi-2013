function varargout = plot_embryo(soln,midline_channel)
%This function plots our embryos.
%
%function varargout = plot_embryo(soln,midline_channel)
%
% "soln": output of "analyze_xs.m" and probably "fit_gaussian.m" too.
% "midline_channel": each channel has the opportunity to have an estimate
%	of the midline.  This index will tell us which channel to use as the
%	estimate for the ventral midline.
%
% "varargout" can be the following three things, in this order:
% - If one output is asked for, varargout = {h_plot}
% - If two are asked for, varargout = {h_plot,leyenda}
% - If three, varargout = {h_plot,leyenda,soln_out}
%
% "h_plot" is an nLSM-by-1 cell array (where nLSM is the number of embryos
%	in "soln") where each element corresponds to all of the object handles
%	of that plot.
% "leyenda": also a cell array similar to "h_plot", but this contains all
%	of the names for each object in each plot.  If you'd like a legend to
%	appear on plot "i", say "legend(leyenda{i})"
% "soln_out": an output version of your input variable "soln". The only
%	thing that might change between "soln" and "soln_out" is if you found
%	the midline by hand.  If that's the case, then "soln_out" will have the
%	newly-found ventral midline values stored in every channel of s_mid.

if ~exist('midline_channel','var')
	midline_channel = NaN;
end

nLSM = length(soln);
h_plot = cell(nLSM,1);
leyenda = h_plot;
if nargout == 3
	soln_out = soln;
end
for i = 1:nLSM
	data = soln(i);
	
	%
	% Determining what to do with the midline.
	%
	alignbyhand = false;
	if isnan(midline_channel) || ...
			midline_channel == 0 || ... 
			isnan(data.s_mid(midline_channel))
		if all(isnan(data.s_mid)) || midline_channel == 0
			midline_channel = 1;
			data.s_mid = zeros(size(data.s_mid));
			alignbyhand = true;
		else
			mc = find(~isnan(data.s_mid));
			midline_channel = mc(1); % default to first non-NaN channel
		end
	end
	
	%
	% Colors
	%
	C = {'b','r',[0 0.8 0],[0.8 0.4 0],[0 0.7 0.7],[0.7 0.7 0]};
	
	figger
% 	set(gcf,'Position',[986   360   672   504])
	set(gcf,'Position',[450  208  672  504])
	fignum = gcf;
	xlim([-1 1])
	box on
	
	[h1,leyenda1,C] = plot_nucprotein(data,fignum,midline_channel,C);
	[h2,leyenda2,C] = plot_mrna(data,fignum,midline_channel,C);
	[h3,leyenda3] = plot_intronic(data,fignum,midline_channel,C);
	
	leyenda{i} = [leyenda1,leyenda2,leyenda3];
	set(gca,'Fontsize',16)
	legend(leyenda{i})
	xlabel('DV coordinate')
	ylabel('Intensity (AU)')
	hold off
	
	h_plot{i} = [h1;h2;h3];
	
	%
	% If there is no known ventral midline, this allows the user to
	% interact with the plot and choose the ventral midline.  A left-click
	% will choose the ventral midline, while any other kind of click or
	% keypress will exit with the current ventral midline chosen.
	%
	if alignbyhand
		s_mid = 0;
		button = 1;
		while button == 1
			[sstar,y,button] = ginput(1);
			if button ~= 1
				break
			end
			
			%
			% Realigning each plot
			%
			for j = 1:length(h_plot{i})
				g = get(h_plot{i}(j));
				
				if isfield(g,'UData') || length(g.YData) ~= length(data.s)
					Xdata = g.XData;
					Xdata = mod(Xdata-sstar+1,2)-1;
					set(h_plot{i}(j),'XData',Xdata)					
				else
					Ydata = g.YData;
					Ydata = circShiftDU(Ydata',sstar)';
					set(h_plot{i}(j),'YData',Ydata)
				end
			end	
			s_mid = s_mid + sstar;
		end	
		s_mid = mod(s_mid+1,2)-1;
		if nargout == 3
			soln_out(i).s_mid = ones(size(data.s_mid))*s_mid;
		end
	end
	
end

switch nargout
	case 1
		varargout = {h_plot};
	case 2
		varargout = {h_plot,leyenda};
	case 3
		varargout = {h_plot,leyenda,soln_out};
	otherwise
		varargout = {};
end

% --------- subfunction to plot the nuclear protein concentrations ------
function [h,leyenda,C] = plot_nucprotein(data,fignum,midline_channel,C)

channels = data.channels;
channels0 = channels;
% channels(channels == 5) = []; % removing the "N/A" channels
channelnames = data.channelnames;
% channelnames(channels0 == 5) = []; % removing the "N/A" channels
j_np = find(channels == 3);
n_np = length(j_np);
gaussnames = data.nucprotein_names;

if n_np ~= 0
	s_mid = data.s_mid(midline_channel);
	s_mid1 = s_mid - data.s_mid(j_np);
	s_mid1 = mod(s_mid1+1,2) - 1;
		
	%
	% Defining s,t
	%
	s = data.S;
	s = mod(s-s_mid+1,2) - 1;
	t = data.R;
	stdev = data.Std_R;
	
	%
	% Plotting each nuclear protein
	%
	figure(fignum)
	hold on
	h1 = errorbar(repmat(s,1,size(t,2)),t,stdev,'.')';
	
	%
	% Changing colors
	%
	for i = 1:n_np
		set(h1(i),'Color',C{i})
	end
	C(1:n_np) = [];
	
	leyenda1 = channelnames(j_np);
	if  isfield(data,'A')  && all(~isnan(data.A))
		leyenda2 = gaussnames;
		
		%
		% Loop through the fit of each nuclear protein
		%
		n_fit = 301;
		s_fit = linspace(-1,1,n_fit)';
		h2 = zeros(n_np,1);
		for j = 1:length(data.A)
			t_fit = data.A(j).*exp(-s_fit.^2/data.sig(j)^2/2) +...
				data.B(j) + data.M(j)*abs(s_fit);
            count = j;
            while ~strcmp(leyenda1{count},gaussnames(j))
                count = count + 1;
            end
			t_fit = circShiftDU(t_fit,s_mid1(count));
			h2(j) = plot(s_fit,t_fit,'k','Linewidth',2);
			leyenda2{j} = [leyenda2{j} ' fit'];
		end
	else
		h2 = [];
		leyenda2 = {};
	end
	
	h = [h1;h2];
	leyenda = [leyenda1 leyenda2];
else
	h = [];
	leyenda = {};
end

% --------- subfunction to plot the mRNA concentrations ------
function [h,leyenda,C] = plot_mrna(data,fignum,midline_channel,C)

channels = data.channels;
channels0 = channels;
channels(channels == 5) = []; % removing the "N/A" channels
channelnames = data.channelnames;
channelnames(channels0 == 5) = [];
jmrna = find(channels == 2);
% jmrna = find(channels == 2 | channels == 3);
if ~isempty(jmrna)
	s_mid = data.s_mid(midline_channel);
	s = data.s;
	t = data.t(:,jmrna);
	
	figure(fignum)
	Ylim = ylim;
	A = 0.9*diff(Ylim);
	B = Ylim(1) + 0.1*diff(Ylim);
% 	if isfield(data.metadata,'genes') && isnumeric(data.metadata.genes.rbr)
% 		rbr = data.metadata.genes.rbr(jmrna);
% 		t = subtrbkgrnd(t,rbr);
% 	end
	t = circShiftDU(t,s_mid);
	hold on
	h = plot(s,A*mDU(t,'col')+B,'Linewidth',2);
	leyenda = channelnames(jmrna);	
	for i = 1:length(h)
		set(h(i),'Color',C{i})
	end
	C(1:length(h)) = [];
else
	h = [];
	leyenda = {};
end


% --------- subfunction to plot the intronic probe concentrations ------
function [h,leyenda] = plot_intronic(data,fignum,midline_channel,C)

channels = data.channels;
channels0 = channels;
channels(channels == 5) = []; % removing the "N/A" channels
channelnames = data.channelnames;
channelnames(channels0 == 5) = [];
j_intron = find(channels == 4);

if ~isempty(j_intron)
	s_mid = data.s_mid(midline_channel);
	s = data.S;
	s = mod(s-s_mid+1,2) - 1;
	t = data.Intron;
	
	figure(fignum)
	Ylim = ylim;
	A = 0.9*diff(Ylim);
	B = Ylim(1) + 0.1*diff(Ylim);
	h = plot(s,A*mDU(t)+B,'.');
	leyenda = channelnames(j_intron);
	for k = 1:length(leyenda)
		leyenda{k} = [leyenda{k} ' intron'];
	end
	for i = 1:length(h)
		set(h(i),'Color',C{i})
	end
else
	h = [];
	leyenda = {};
end



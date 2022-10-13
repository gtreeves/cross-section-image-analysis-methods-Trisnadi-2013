function soln_out = fit_gaussian(soln)
%Fits the data of the dl nuclear gradient to a gaussian.
%
%function soln_out = fit_gaussian(soln)
%
% This function accepts as input the structure "soln", which is the output
% of the function "analyze_xs".  That structure contains intensity values
% of the nuclei and of nuclear protein (e.g., Dorsal).  This function's
% main purpose is to approximate the intensity profile of nuclear proteins
% as gaussians.
%
% Inputs:
%
% "soln": structure; output from analyze_xs
%
% Outputs:
%
% "soln_out": a structure identical to the input, "soln", but also contains
% the following extra fields:  
%	"A": gradient amplitude
%	"B": basal levels
%	"M": the outer slope
%	"mu": the midline associated with the nuclear protein channel
%	"sig": the spatial extent of the gradient
%	"dA,dB,dM,dmu,dsig": errorbars on these estimates
%	"gof": the R-squared goodness of fit for each nuclear protein channel
%
% Also, some more metadata is included under the soln.metadata field.

warning off all
soln_out = soln;
nLSM = length(soln_out);

%
% Going through each embryo.
%
for i = 1:nLSM
	data = soln(i);
	channels = data.channels;
	channelnames = data.channelnames;
	channels0 = channels;
% 	channels(channels0 == 5) = []; % removing the "N/A" channels
% 	channelnames(channels0 == 5) = [];
	
	k_np = find(channels == 3); % indices of channels with nuc proteins
	n_np = length(k_np);
	
	%
	% Finding if dl (or any other gaussian-like protein) is present in
	% "channelnames".  If so, run "fit_gaussian".  If not, move on.  It is
	% possible to put multiple nuclear proteins in the same channel.  in
	% that case, all proteins in the channel must be gaussian-like for this
	% to work.  Possible example: dl and pMad in the same channel.  The
	% line below, which defines the variable "gaussianmatches", contains a
	% list of all nuclear proteins whose distribution is expected to match
	% the gaussian profile.  Edit this list if you wish more or fewer
	% proteins to be sent to "fit_gaussian".
	%
	gaussianmatches = {'dl','pMad','dl-GFP','dl-Venus','Med-GFP'};
	cnames = channelnames(k_np);	
	for j = 1:n_np % loop through all nuc protein channels
		cname = cnames{j};
		cname2 = str2cell(cname,',')'; % split channelname in case mult
		% proteins present.
		
		v = false(length(cname2),1);
		for k = 1:length(cname2)
			v(k) = ~isempty(strmatch(cname2{k},gaussianmatches,'exact'));
		end
		if ~all(v)
			k_np(j) = NaN;
		end
		
	end
	
	k_np0 = k_np; n_np0 = n_np; % includes non-gaussian nucproteins
	k_np(isnan(k_np)) = []; % removing channels that aren't all gaussian	
	n_np = length(k_np);
	
	if n_np > 0
		A = NaN(1,n_np);
		B = A; M = A; sig = A; mu = A; R2 = A;
		dA = A; dB = A; dM = A; dmu = A; dsig = A;
		cint = cell(1,n_np); cint68 = cint;
		
		%
		% Put together vectors s,r
		%
		s = data.S;
		r = data.R;
		stdev = data.Std_R;
		
		%
		% Going through each nuclear protein
		%
		for j = 1:n_np
			
			r1 = r(:,k_np0 == k_np(j));
			stdev1 = stdev(:,k_np0 == k_np(j));
			[A1,B1,M1,mu1,sig1,R21,cint1,cint681] = fitelephant(s,r1,stdev1);
			A(j) = A1;
			B(j) = B1;
			M(j) = M1;
			mu(j) = mu1;
			sig(j) = sig1;
			
			dc = diff(cint681)/2;
			dA(j) = dc(1);
			dB(j) = dc(2);
			dM(j) = dc(3);
			dmu(j) = dc(4);
			dsig(j) = dc(5);
			
			R2(j) = R21;
			cint{j} = cint1;
			cint68{j} = cint681;
		end
		
		soln_out(i).nucprotein_names = channelnames(k_np);
		
		soln_out(i).A = A;
		soln_out(i).B = B;
		soln_out(i).M = M;
		soln_out(i).mu = mu;
		soln_out(i).sig = sig;
		
		soln_out(i).dA = dA;
		soln_out(i).dB = dB;
		soln_out(i).dM = dM;
		soln_out(i).dmu = dmu;
		soln_out(i).dsig = dsig;
		
		soln_out(i).gof = R2;
		
		soln_out(i).metadata.cint = cint;
		soln_out(i).metadata.cint68 = cint68;
		
		soln_out(i).s_mid(k_np) = mean(mu,1);
		% there is an s_mid for each channel.
	end
end

% --------------- subfunction to perform the fit ---------------
function [A,B,M,mu,sig,R2,cint,cint68] = fitelephant(s,t,stdev)

%
% Getting an estimate of where the max is
%
[ssort,isort] = sort(s); % sorting wrt s
tsort = t(isort);
p = 10; % periodic extension
t_ext = [tsort(end-p+1:end);tsort;tsort(1:p)];
s_ext = [ssort(end-p+1:end)-2;ssort;ssort(1:p)+2];
tsmooth = smooth(s_ext,t_ext,p); % smoothing
tsmooth = tsmooth(p+1:end-p); 
[tmax,imax] = max(tsmooth); tmin = min(tsmooth);

%
% Shifting so that the max of the dl grad is near the middle;
% that is, near "s = 0;".
%
s_mid1 = ssort(imax);
s1 = mod(s-s_mid1+1,2)-1;

%
% The initial guesses and upper and lower bounds of each parameter
%
A = tmax - tmin; AL = 0.1*A; AU = 10*A;
B = tmin; BL = 0; BU = A + tmin;
M = 0; ML = -1e6; MU = 1e6;
mu = 0; muL = -0.3; muU = 0.3;
sig = 0.15; sigL = 0.01; sigU = 1;

%
% Defining options and the equation to fit to
%
f = fittype('A*exp(-(x-mu)^2/2/sigma^2)+B+M*abs(x-mu)');
opts = fitoptions('Method','NonlinearLeastSquares',...
	'Startpoint',[A B M mu sig],...
	'Lower',[AL BL ML muL sigL],...
	'Upper',[AU BU MU muU sigU],...
	'Weights',1./stdev);

%
% The actual fit
%
[cfun,gof] = fit(s1,t,f,opts);
R2 = gof.rsquare;

%
% Parameter values
%
coeffvals = coeffvalues(cfun);
A = coeffvals(1);
B = coeffvals(2);
M = coeffvals(3);
s_mid2 = coeffvals(4);
sig = coeffvals(5);
mu = s_mid1 + s_mid2;

%
% Confidence intervals for each parameter
%
cint = confint(cfun);
cint68 = confint(cfun,(1-2*(1-normcdf(1,0,1))));
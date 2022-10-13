function [Y,stdY] = intronicintensity(I,nucstats,intr_pad)
%Computes intronic probe intensity of given channels.  Works w/analyze_xs.
%
%function [Y,stdY] = intronicintensity(I,nucstats)
%
% This function calculates the intronic probe (nuclear dot) intensity for
% each nucleus in I, specified by the pixel index lists in "nucstats".  The
% structure "nucstats" comes from a call to "regionprops" from
% "find_nuclei". 
%
% The nuclear dot intensity is taken to be the median intensity of an
% n-by-n neighborhood centered around the max intensity pixel (which is
% presumably where the nuclear dot actually lies), for each nucleus.  The
% reason why we take the median of a neighborhood is to decrease the noise
% we get from randomly-activated pixels (given most of these images are
% derived from PMT's).
%
% The value of "n" is "2*intr_pad+1", where "intr_pad" is equal to 2
% (although this could be adjusted directly in the code below if you find
% this is too much of a padding).
%
% The image "I" is a z-slice of an embryo cross
% section, and only contains the channels that have intronic probe.  The
% order is specified by the original image; so channels are removed from I
% but the remaining ones are not shuffled around.
%
% Inputs:
%
% "I": image of z-slice of embryo cross section, including only the
%	channels that have intronic probe.
% "nucstats": structure from "find_nuclei" (calculated by "regionprops")
%	that contains the PixelList and PixelIdxList of each nucleus.
% "intr_pad": the padding for median filter on the intronic dots.  This is
%	an optional argument, and its default value (if not passed) is 2.
%
% Outputs:
% 
% "Y": n_nuc-by-numch array of mean nuclear intensities.
% "stdY": misnomer, should be semY since it's the standard error of the
%	mean nuclear intensity.

if ~exist('intr_pad','var')
	intr_pad = 2; % padding for median filter on the intronic dots.
end
[H,W,numch] = size(I);
n_nuc = length(nucstats);
Y = zeros(n_nuc,numch); stdY = Y;

%
% Run a loop to go through each nucleus
%
for i = 1:length(nucstats)
	v = nucstats(i).PixelIdxList; % the pixels on nucleus "i"
	if ~isempty(v)
		for j = 1:numch			
			I1 = I(:,:,j);
			[Imax,imax] = max(I1(v)); % find the max intensity pixel
			
			%
			% Converting the pixel index into an x-pixel,y-pixel pair.
			%
			px = nucstats(i).PixelList(imax,:);
			ii = px(2); jj = px(1); % PixelList gives x first, then y.
			
			%
			% In case the max intensity pixel is closer than "intr_pad"
			% distance from the edge of the image
			%
			i1 = max(ii-intr_pad,1); i2 = min(ii+intr_pad,H);
			j1 = max(jj-intr_pad,1); j2 = min(jj+intr_pad,W);
			
			%
			% The n-by-n neighborhood
			%
			J = double(I1(i1:i2,j1:j2));
			
			Y(i,j) = median(J(:));
			stdY(i,j) = (max(J(:)) - min(J(:)))/2;
		end
	end
end
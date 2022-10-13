function s = lsmReadDU(varargin)
%Ftn to read in image z-stacks from lsm fmt, from tiffread ver 2.5
%
%function s = lsmReadDU(varargin)
%
% Reads 8,16,32 bits uncompressed grayscale and color tiff stacks out of
% the lsm format from Zeiss.  The entire TIFF standard is not supported,
% but you may extend it by editing this file. If you do so, it would be
% nice to return your modification to F. Nedelec, so that it can be
% included in future releases.
%
% The function can be called with a file name in the current directory, or
% without argument, in which case it pops up a file openning dialog to
% allow manual selection of the file. If the stack contains multiple
% images, loading can be restricted by specifying the first and last images
% to read, or just one image to read.
%
% At return, "img_read" contains the number of images read, and "s" is a
% structure containing the different images with some additional
% information. The images themselves are stored in the field .data, a cell
% array.  Each element of the cell array is a z-slice image, and is stored
% as a 3D array.  Note that the identity of RGB may be switched around. The
% pixel values are returned in the native (integer) format, and must be
% converted to be used in most matlab functions.
%
% Example:
% im = tiffread('spindle.lsm');
% imshow(im.data{1});
%
% Optional argument varargin can consist of the following things:
%	(1) "filename": string containing the filename (including relative
%		path) to the lsm file.  If not specified, a dialog box will appear
%		for user to browse to the file.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	(2,3) "index1,index2": If you wish to extract only a specific z-slice,
%		specify only "index1".  If you wish to extract an inverval of
%		slices, specify the fist in "index1", and the second in "index2".
%
% "s": structure containing the information about the images, as well as
% the images themselves in the field ".data": 
%	"data": field of "stack", which is a cell array.  Each element of the
%		cell array is a z-slice of the image in the lsm file.
%
% Greg Reeves, Caltech, Feb 2008.  Adapted to read in the whole z-stack of
% lsm files, from:
%
% Francois Nedelec, EMBL, Copyright 1999-2007.
% rewriten July 7th, 2004 at Woods Hole during the physiology course.
% last modified December 21, 2007.
% With contributions from:
%   Kendra Burbank for the waitbar
%   Hidenao Iwai for the code to read floating point images,
%   Stephen Lang to be more compliant with PlanarConfiguration
%   Jan-Ulrich Kreft for Zeiss LSM support
%   Elias Beauchanp and David Kolin for additional Metamorph support
%
% Please, help us improve this software: send us feedback/bugs/suggestions
% This software is provided at no cost by a public research institution.
%
% Francois Nedelec
% nedelec (at) embl.de
% Cell Biology and Biophysics, EMBL; Meyerhofstrasse 1; 69117 Heidelberg;
%	Germany
% http://www.embl.org
% http://www.cytosim.org

%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	filename = varargin{iArg}; else
	% if this isn't specified, the user browses to the file:
	[filename, pathname] = ...
		uigetfile('*.tif;*.stk;*.lsm', 'select image file');
	filename = [ pathname, filename ];
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	index1 = varargin{iArg}; else
	%
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	index2 = varargin{iArg}; else
	%
end, iArg = iArg + 1;

if exist('index1','var') && ~exist('index2','var')
	index2 = index1;
elseif ~exist('index1','var') && ~exist('index2','var')
	index1 = 1; index2 = Inf;
end

%
% Preamble.  The structure IMG is returned to the user, while TIF is not.
% so tags usefull to the user should be stored as fields in IMG, while
% those used only internally can be stored in TIF.
%
TIF.SampleFormat = 1; % setting default values
TIF.SamplesPerPixel = 1;
TIF.BOS = 'ieee-le'; % byte order string

TIF.file = fopen(filename,'r','l');
if TIF.file == -1
	error(['file <',filename,'> not found.']);
end

%
% Reading in header, including byte order (II = little endian, MM = big
% endian), whether we have tif images (id must be 42), and then the byte
% offset for the first image file directory (IFD)
%
byte_order = fread(TIF.file, 2, '*char');
if ( strcmp(byte_order', 'II') )
	TIF.BOS = 'ieee-le'; % normal PC format
elseif ( strcmp(byte_order','MM') )
	TIF.BOS = 'ieee-be';
else
	error('This is not a TIFF file (no MM or II).');
end

tiff_id = fread(TIF.file,1,'uint16', TIF.BOS);
if (tiff_id ~= 42) % if it's 42, then it's TIFF fmt.
	error('This is not a TIFF file (missing 42).');
end

TIF.img_pos = fread(TIF.file, 1, 'uint32', TIF.BOS); % byte offset, 1st IFD

nImages = 1; % the nth image we're on.
while TIF.img_pos ~= 0
	clear IMG
	IMG.filename = fullfile(pwd,filename);

	%
	% Reading in all entries in this IFD.
	%
	[TIF,IMG] = readIFDentries(TIF,IMG);

% 	%
% 	% Total number of bytes per image:
% 	%
% 	PlaneBytesCnt = IMG.width*IMG.height*TIF.BytesPerSample;
	
	%
	% Reading in the next IFD address:
	%
    TIF.img_pos = fread(TIF.file, 1, 'uint32', TIF.BOS);

	n_slices = IMG.lsm.DimensionZ;

	for i = index1:min(n_slices,index2)

		TIF.StripCnt = 1;
		if TIF.SamplesPerPixel == 1
			planeCnt = 1; [TIF,IMG.data{i}] = read_plane(TIF,IMG,planeCnt);
		else
			planeCnt=1;
			[TIF,IMG.data{i}(:,:,1)] = read_plane(TIF,IMG,planeCnt);
			planeCnt=2;
			[TIF,IMG.data{i}(:,:,2)] = read_plane(TIF,IMG,planeCnt);
			planeCnt=3;
			[TIF,IMG.data{i}(:,:,3)] = read_plane(TIF,IMG,planeCnt);
		end

	end
	nImages = nImages + 1;
end

%
% Cleaning up
%
fclose(TIF.file);
if exist('waitbar_handle', 'var')
	delete( waitbar_handle );
	clear waitbar_handle;
end
s = IMG;


% --------- subfunction to read the image plane (z-slice) --------
function [TIF,plane] = read_plane(TIF,IMG,planeCnt)

offset = 0;
width = IMG.width;
height = IMG.height;

%
% Returning an empty array if the sample format has zero bits
%
if TIF.BitsPerSample(planeCnt) == 0
	plane = [];
	return;
end

% fprintf(1,'reading plane %i size %i %i\n',planeCnt,width,height);

%
% Determining the type needed to store the pixel values:
%
switch TIF.SampleFormat
	case 1
		classname = sprintf('uint%i', TIF.BitsPerSample(planeCnt));
		% This would be uint8, I think.
	case 2
		classname = sprintf('int%i', TIF.BitsPerSample(planeCnt));
	case 3
		if TIF.BitsPerSample(planeCnt) == 32
			classname = 'single';
		else
			classname = 'double';
		end
	otherwise
		error('unsuported TIFF sample format %i', TIF.SampleFormat);
	%
end

%
% Preallocate a matrix to hold the sample data:
%
plane = zeros(width,height,classname);

%
% Reading the strips and concatenating them:
%
line = 1;
while TIF.StripCnt <= TIF.StripNumber

	[TIF,strip] = read_strip(TIF,offset,width,planeCnt,classname);
	TIF.StripCnt = TIF.StripCnt + 1;

	%
	% Copying the strip onto the data
	%
	plane(:,line:(line+size(strip,2)-1)) = strip;

	line = line + size(strip,2);
	if line > height
		break
	end

end

%
% Extract valid part of data if needed
%
if ~all(size(plane) == [width height])
	plane = plane(1:width,1:height);
	warning('Cropping data: found more bytes than needed...')
end

%
% transposing the image (otherwise display is rotated in matlab)
%
plane = plane';



% --------- subfunction to read strip --------
function [TIF,strip] = ...
	read_strip(TIF,offset,width,planeCnt,classname)

% fprintf(1,'reading strip at position %i\n',...
%	TIF.StripOffsets(stripCnt)+offset)
stripCnt = TIF.StripCnt
StripLength = TIF.StripByteCounts(stripCnt)./TIF.BytesPerSample(planeCnt);

% fprintf(1,'reading strip %i\n',stripCnt);
fseek(TIF.file,TIF.StripOffsets(stripCnt)+offset,-1);
bytes = fread(TIF.file,StripLength,classname,TIF.BOS );

if any(length(bytes) ~= StripLength)
	error('End of file reached unexpectedly.');
end

strip = reshape(bytes,width,StripLength/width);




% --------- subfunction to store the type of data --------
function [nbBytes, matlabType] = convertType(tiffType)

switch (tiffType)
	case 1
		nbBytes=1;
		matlabType='uint8';
	case 2
		nbBytes=1;
		matlabType='uchar';
	case 3
		nbBytes=2;
		matlabType='uint16';
	case 4
		nbBytes=4;
		matlabType='uint32';
	case 5
		nbBytes=8;
		matlabType='uint32';
	case 11
		nbBytes=4;
		matlabType='float32';
	case 12
		nbBytes=8;
		matlabType='float64';
	otherwise
		error('tiff type %i not supported', tiffType)
	%
end




% --------- subfunction to read IFD entry --------
function  [TIF,IMG] = readIFDentries(TIF,IMG)

%
% Move in the file to the beginning of the IFD
%
fseek(TIF.file, TIF.img_pos, -1);
% disp(strcat('reading img at pos :',num2str(TIF.img_pos)));

%
% Reading in the number of IFD entries
%
num_entries = fread(TIF.file,1,'uint16', TIF.BOS);
% disp(strcat('num_entries =', num2str(num_entries)));

%
% Read and process each IFD entry
%
for i = 1:num_entries
	
	file_pos  = ftell(TIF.file); % saving current position
	TIF.entry_tag = ...
		fread(TIF.file,1,'uint16',TIF.BOS); % read entry tag

	entry.tiffType = fread(TIF.file,1,'uint16',TIF.BOS);
	entry.cnt = fread(TIF.file,1,'uint32',TIF.BOS);
	% disp(['tiffType =',num2str(entry.tiffType),...
	%	',cnt = ',num2str(entry.cnt)])

	[entry.nbBytes,entry.matlabType] = convertType(entry.tiffType);

	if entry.nbBytes*entry.cnt > 4
		% next field contains an offset:
		offset = fread(TIF.file,1,'uint32',TIF.BOS);
		% disp(strcat('offset = ',num2str(offset)));
		fseek(TIF.file,offset,-1);
	end

	if TIF.entry_tag == 34412
		entry.val = readLSMinfo(TIF);
	else
		if entry.tiffType == 5
			entry.val=fread(TIF.file,2*entry.cnt,entry.matlabType,TIF.BOS);
		else
			entry.val = fread(TIF.file,entry.cnt,entry.matlabType,TIF.BOS);
		end
	end

	if entry.tiffType == 2;
		entry.val = char(entry.val');
	end

	% disp(strcat('reading entry <',num2str(TIF.entry_tag),'>'));

	%
	% Depending on the value of the entry tag (stored in TIF.entry_tag),
	% the structure "entry", which was read in from the lsm file using the
	% subfunction "readIFDentry", will mean different things.  Here, we
	% unpack those different meanings.
	%
	switch TIF.entry_tag
		case 254
			TIF.NewSubfiletype = entry.val;
		case 256 % image width - number of columns
			IMG.width = entry.val;
		case 257 % image height - number of rows
			IMG.height = entry.val;
			TIF.ImageLength = entry.val;
		case 258 % BitsPerSample per sample
			TIF.BitsPerSample = entry.val;
			TIF.BytesPerSample = TIF.BitsPerSample / 8;
			IMG.bits = TIF.BitsPerSample(1);
			% fprintf(1,'BitsPerSample %i %i %i\n', entry.val);
		case 259 % compression
			if entry.val ~= 1
				error('Compression format not supported.');
			end
		case 262 % photometric interpretation
			TIF.PhotometricInterpretation = entry.val;
			if ( TIF.PhotometricInterpretation == 3 )
				warning('Ignoring TIFF look-up table');
			end
		case 269
			IMG.document_name = entry.val;
		case 270 % comments:
			IMG.info = entry.val;
		case 271
			IMG.make = entry.val;
		case 273 % strip offset
			TIF.StripOffsets = entry.val;
			TIF.StripNumber = entry.cnt;
			% fprintf(1,'StripNumber = %i,...
			%	size(StripOffsets) = %i %i\n',...
			%	TIF.StripNumber,size(TIF.StripOffsets));
		case 277 % sample_per pixel
			TIF.SamplesPerPixel = entry.val;
			% fprintf(1,'Color image: sample_per_pixel=%i\n',...
			%	TIF.SamplesPerPixel);
		case 278 % rows per strip
			TIF.RowsPerStrip = entry.val;
		case 279 % strip byte counts - number of bytes in each strip
			% after any compression
			TIF.StripByteCounts = entry.val;
		case 282 % X resolution
			IMG.x_resolution = entry.val;
		case 283 % Y resolution
			IMG.y_resolution = entry.val;
		case 284 % planar configuration describe the order of RGB
			TIF.PlanarConfiguration = entry.val;
		case 296 % resolution unit
			IMG.resolution_unit = entry.val;
		case 305 % software
			IMG.software = entry.val;
		case 306 % datetime
			IMG.datetime = entry.val;
		case 315
			IMG.artist = entry.val;
		case 317 % predictor for compression
			if entry.val ~= 1
				error('unsuported predictor value');
			end
		case 320 % color map
			IMG.cmap = entry.val;
			IMG.colors = entry.cnt/3;
		case 339
			TIF.SampleFormat = entry.val;
		case 34412 % Zeiss LSM data
			IMG.lsm = entry.val;
		otherwise
			fprintf(1,'Ignored TIFF entry with tag %i (cnt %i)\n',...
				TIF.entry_tag,entry.cnt);
			%
	end


	fseek(TIF.file, file_pos+12, -1); % move to next IFD entry
end


% --------- subfunction to parse LSM info --------
function R = readLSMinfo(TIF)

%
% Only the first table is read! This provides only very partial info, as
% the offset indicates that additional data is stored in the file.
%

R.MagicNumber = sprintf('0x%x',fread(TIF.file, 1, 'uint32', TIF.BOS));
StructureSize = fread(TIF.file, 1, 'int32', TIF.BOS);
R.DimensionX = fread(TIF.file, 1, 'int32', TIF.BOS);
R.DimensionY = fread(TIF.file, 1, 'int32', TIF.BOS);
R.DimensionZ = fread(TIF.file, 1, 'int32', TIF.BOS);
R.DimensionChannels = fread(TIF.file, 1, 'int32', TIF.BOS);
R.DimensionTime = fread(TIF.file, 1, 'int32', TIF.BOS);
R.DataType = fread(TIF.file, 1, 'int32', TIF.BOS);
R.ThumbnailX = fread(TIF.file, 1, 'int32', TIF.BOS);
R.ThumbnailY = fread(TIF.file, 1, 'int32', TIF.BOS);
R.VoxelSizeX = fread(TIF.file, 1, 'float64', TIF.BOS);
R.VoxelSizeY = fread(TIF.file, 1, 'float64', TIF.BOS);
R.VoxelSizeZ = fread(TIF.file, 1, 'float64', TIF.BOS);
R.ScanType = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DataType = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetVectorOverlay = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetInputLut = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetOutputLut = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetChannelColors = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.TimeInterval = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OffsetChannelDataTypes = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetScanInformation = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetKsData = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetTimeStamps = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetEventList = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetRoi = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetBleachRoi = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.OffsetNextRecording = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.Reserved = fread(TIF.file, 1, 'uint32', TIF.BOS);

% 
% % --------- subfunction to read IFD entry --------
% function  [TIF,IMG] = readIFDentry(TIF,IMG)
% 
% entry.tiffType = fread(TIF.file, 1, 'uint16', TIF.BOS);
% entry.cnt      = fread(TIF.file, 1, 'uint32', TIF.BOS);
% % disp(['tiffType =',num2str(entry.tiffType),',cnt = ',num2str(entry.cnt)])
% 
% [entry.nbBytes,entry.matlabType] = convertType(entry.tiffType);
% 
% if entry.nbBytes * entry.cnt > 4
% 	% next field contains an offset:
% 	offset = fread(TIF.file,1,'uint32',TIF.BOS);
% 	% disp(strcat('offset = ',num2str(offset)));
% 	fseek(TIF.file,offset,-1);
% end
% 
% if TIF.entry_tag == 34412
% 	entry.val = readLSMinfo(TIF);
% else
% 	if entry.tiffType == 5
% 		entry.val = fread(TIF.file,2*entry.cnt,entry.matlabType,TIF.BOS);
% 	else
% 		entry.val = fread(TIF.file,entry.cnt,entry.matlabType,TIF.BOS);
% 	end
% end
% 
% if entry.tiffType == 2;
% 	entry.val = char(entry.val');
% end
% 
% % disp(strcat('reading entry <',num2str(TIF.entry_tag),'>'));
% 
% %
% % Depending on the value of the entry tag (stored in TIF.entry_tag),
% % the structure "entry", which was read in from the lsm file using the
% % subfunction "readIFDentry", will mean different things.  Here, we
% % unpack those different meanings.
% %
% switch TIF.entry_tag
% 	case 254
% 		TIF.NewSubfiletype = entry.val;
% 	case 256 % image width - number of columns
% 		IMG.width = entry.val;
% 	case 257 % image height - number of rows
% 		IMG.height = entry.val;
% 		TIF.ImageLength = entry.val;
% 	case 258 % BitsPerSample per sample
% 		TIF.BitsPerSample = entry.val;
% 		TIF.BytesPerSample = TIF.BitsPerSample / 8;
% 		IMG.bits = TIF.BitsPerSample(1);
% 		% fprintf(1,'BitsPerSample %i %i %i\n', entry.val);
% 	case 259 % compression
% 		if entry.val ~= 1
% 			error('Compression format not supported.');
% 		end
% 	case 262 % photometric interpretation
% 		TIF.PhotometricInterpretation = entry.val;
% 		if ( TIF.PhotometricInterpretation == 3 )
% 			warning('Ignoring TIFF look-up table');
% 		end
% 	case 269
% 		IMG.document_name = entry.val;
% 	case 270 % comments:
% 		IMG.info = entry.val;
% 	case 271
% 		IMG.make = entry.val;
% 	case 273 % strip offset
% 		TIF.StripOffsets = entry.val;
% 		TIF.StripNumber = entry.cnt;
% 		% fprintf(1,'StripNumber = %i,...
% 		%	size(StripOffsets) = %i %i\n',...
% 		%	TIF.StripNumber,size(TIF.StripOffsets));
% 	case 277 % sample_per pixel
% 		TIF.SamplesPerPixel = entry.val;
% 		% fprintf(1,'Color image: sample_per_pixel=%i\n',...
% 		%	TIF.SamplesPerPixel);
% 	case 278 % rows per strip
% 		TIF.RowsPerStrip = entry.val;
% 	case 279 % strip byte counts - number of bytes in each strip
% 		% after any compression
% 		TIF.StripByteCounts = entry.val;
% 	case 282 % X resolution
% 		IMG.x_resolution = entry.val;
% 	case 283 % Y resolution
% 		IMG.y_resolution = entry.val;
% 	case 284 % planar configuration describe the order of RGB
% 		TIF.PlanarConfiguration = entry.val;
% 	case 296 % resolution unit
% 		IMG.resolution_unit = entry.val;
% 	case 305 % software
% 		IMG.software = entry.val;
% 	case 306 % datetime
% 		IMG.datetime = entry.val;
% 	case 315
% 		IMG.artist = entry.val;
% 	case 317 % predictor for compression
% 		if entry.val ~= 1
% 			error('unsuported predictor value');
% 		end
% 	case 320 % color map
% 		IMG.cmap = entry.val;
% 		IMG.colors = entry.cnt/3;
% 	case 339
% 		TIF.SampleFormat = entry.val;
% 	case 34412 % Zeiss LSM data
% 		IMG.lsm = entry.val;
% 	otherwise
% 		fprintf(1,'Ignored TIFF entry with tag %i (cnt %i)\n',...
% 			TIF.entry_tag,entry.cnt);
% 		%
% end
% 
% 

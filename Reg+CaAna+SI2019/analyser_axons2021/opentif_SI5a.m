function [header,Aout, imgInfo] = opentif_SI5a(varargin)
%% function [header,Aout, imgInfo] = opentif(varargin)
% Opens a ScanImage TIF file, extracting its header information and, if specified, stores all of image contents as output array Aout if specified. 
% By default, Aout, if specified for output, is of size MxNxCxFxSxV, where C spans the channel indices, F spans the frame indicies, S spans the 
% slice indices, and V the volume indices.
%
% NOTE: IF the second output argument (Aout) is not assigned to output variable
%       THEN image file is not actually read -- only  header information is extracted
%       
% IMPORTANT: opentif currently only exports the header and sequential image data. Once the tiff header specification reaches a stable 
%			 point, parsing and data organization will be reincorporated ++++
%
%% SYNTAX
%   opentif()
%   opentif(filename)
%   header = opentif(...)
%   [header,Aout] = opentif(...)
%   [header,Aout,imgInfo] = opentif(...)
%		INPUT
%       	filename: Name of TIF file, with or without '.tif' extension. If omitted, a dialog is launched to allow interactive selection.
%       	flagN/flagNArg: Flags (string-valued) and/or flag/value pairs, in any order, specifying options to use in opening specified file
%
%		OUTPUT
%       	header: Structure comprising information stored by ScanImage into TIF header
%       	Aout: MxNxCxFxSxV array, with images of size MxN for C channels, F frames, S slices, and V volumes. Default type is uint16. 
%       	imgInfo: Structure comprising basic information about the structure of the output array Aout
%
% NOTE: IF the second output argument (Aout) is not assigned to output variable
%       THEN image file is not actually read -- only header information is extracted
%
%% FLAGS (case-insensitive)
%
%   WITH ARGUMENTS
%       'channel' or 'channels': Argument specifies subset of channel(s) to extract. Ex: 1,[1 3], 2:4. 
%       'frame' or 'frames': Argument specifies subset of frames present to extract. Use 'inf' to specify all frames above highest specified value. Ex: 1:30, [50 inf], [1:9 11:19 21 inf]
%       'slice' or 'slices': Argument specifies subset of slices present to extract. Use 'inf' to specify all slices above highest specified value. Ex: 1:30, [50 inf], [1:9 11:19 21 inf]
%       'volume' or 'volumes': Argument specifies subset of volumes present to extract. Use 'inf' to specify all slices above highest specified value. Ex: 1:30, [50 inf], [1:9 11:19 21 inf]
%
%% NOTES
%   This function replaces the scim_openTif() function supplied with ScanImage 4.2
%  	
%	In case of errors, the program will attempt to output whatever image data is available to it as an uncategorized stream of images
%	This stream will be an array of the form MxNxImg raw ouput without any post-processing, containing all the frames found within the file, where Img is the number of images
%
%   TODO: Port more advanced features to ScanImage 5 from SI3/4 scim_openTif
%   TODO: Add a flag to discard fastZ-flyback frames if present
%

  %% Initialize output variables
  header = [];
  Aout   = [];
  imgInfo = struct();

  %% Constants/Inits
  if nargout < 0 || nargout > 3
      most.idioms.warn('Invalid arguments'); 
      return
  end

  %% Parse input arguments

  flagNames = {'channel' 'channels' 'slice' 'slices' 'frame' 'frames' 'volume' 'volumes'};
  argFlags = {'channel' 'channels' 'slice' 'slices' 'frame' 'frames' 'volume' 'volumes'};

  flagIndices = find(cellfun(@(x)ischar(x) && (ismember(lower(x),flagNames) || ismember(lower(x),argFlags)),varargin));

  flags = cellfun(@lower,varargin(flagIndices),'UniformOutput',false);
  if isempty(flags)
      flags = {};
  end

  streamOutput = false;

  %% Determine input file
  if isempty(find(flagIndices==1)) && nargin>=1 && ischar(varargin{1})
      fileName = varargin{1};
  else
      fileName = '';
  end

  if isempty(fileName)
      [f, p] = uigetfile({'*.tif;*.tiff'},'Select Image File');
      if f == 0
          most.idioms.warn('Invalid arguments'); 
          return;
      end
      fileName = fullfile(p,f); 
  end

  %Extract filepath for future use
  %[filePath,fileStem,fileExt] = fileparts((fileName));

  %% Read TIFF file; extract # frames & image header
  if ~exist(fileName,'file') && ~exist([fileName '.tif'],'file') && ~exist([fileName '.tiff'],'file') 
      error('''%s'' is not a recognized flag or filename. Aborting.',fileName);
  elseif exist([fileName '.tif'],'file') 
      fileName = [fileName '.tif'];
  elseif exist([fileName '.tiff'],'file') 
      fileName = [fileName '.tiff'];
  end

  %most.idioms.warn(['Loading file ' fileName]);

  warning('off','MATLAB:tifflib:TIFFReadDirectory:libraryWarning');
  hTif = Tiff(fileName);

  [header, si_ver, numhdrs, headerLeadStr, numImages] = parseHeaderToStruct(hTif);        

  %Reincorporate conditional once header spec is stable
  if numImages == 0 || isempty(fieldnames(header)) || numhdrs == 0 || si_ver == 0
    quitFlag = streamOutputQuit();
    if quitFlag
      return
    end
  end


  hdr = extractHeaderData(header,si_ver);

  %% Read image meta-data
  savedChans = hdr.savedChans;

  %Display channel information to user
  %most.idioms.warn(['Matrix of channels saved: ' mat2str(savedChans)]);

  numChans = length(savedChans);
  numPixels = hdr.numPixels;
  numLines = hdr.numLines;
  numSlices = hdr.numSlices;
  numVolumes = hdr.numVolumes;
  numFrames = hdr.numFrames;
  numDiscardFrames = 0;
  discardFlybackframesEnabled = false;

  % If using FastZ, use slices value that contains potential flyback frames
  % for proper organization of output image-array
  if si_ver == 2015 && hdr.discardFlybackframesEnabled
    discardFlybackframesEnabled = hdr.discardFlybackframesEnabled;
    % We list the total number of slices so our next trigger checks don't run into trouble
    numSlices = hdr.numFramesPerVolume;   
    % After nextTrigger checking, we'll remove the number of discarded frames from numSlices
    numDiscardFrames = hdr.numDiscardFrames;
  end

  if numSlices > 1 && numFrames > 1
    most.idioms.warn('Cannot interpret multiple frames and slices simultaneously at this time.');
    quitFlag = streamOutputQuit();
    if quitFlag
      return
    end
  end

  %% This section makes sure there are no issues with nextTrigger data
  if numImages ~= numChans*numFrames*numSlices*numVolumes
    % We are working under the assumption that only volumes can have multiple "slices"
    if numSlices > 1
      numVolumes = floor(numImages/numChans/numSlices);
      numFrames = 1;  % This should already be the case
    elseif numFrames > 1
      % In this case there are no volumes, since we only can have 1 frame and multiple slices in a volume
      numVolumes = 1; % This should already be the case
      numSlices = 1;  % This should already be the case
      % We discard the previous value of frames and adjust to what was acquired before the next-trigger came in
      numFrames = floor(numImages/numChans);  
    end

    if numImages ~= numChans*numFrames*numSlices*numVolumes
      most.idioms.warn('Unexpected number of images.');
      quitFlag = streamOutputQuit();
      if quitFlag
        return
      end
    end
  end

  % This takes care of any flybackFrames in the current volume
  if si_ver == 2015
    numSlices = numSlices - numDiscardFrames;
  end

  %DEBUG msg
  %most.idioms.warn(['numImages = ' num2str(numImages)]);
  %most.idioms.warn(['numChans = ' num2str(numChans)]);
  %most.idioms.warn(['numFrames = ' num2str(numFrames)]);
  %most.idioms.warn(['numSlices = ' num2str(numSlices)]);
  %most.idioms.warn(['numVolumes = ' num2str(numVolumes)]);
  %most.idioms.warn(' ');

  if ~numFrames || ~numSlices
    most.idioms.warn('Acquisition did not complete a single frame or slice. Aborting.');
    quitFlag = streamOutputQuit();
    if quitFlag
      return
    end
  end

  %Remove extra headers for multiple channels since they are identical
  numUniqueHeaders = floor(numImages/numChans);  % This should account for ext trigger mode
  if numChans > 1
      tempHeader = zeros(numhdrs, numUniqueHeaders); 
      
      for iter = 1:numhdrs    
          headSame = zeros(numChans, numUniqueHeaders); 

          for i = 1:numChans
              eval(['headSame(i, :) = header.' headerLeadStr{iter} '(1, i:numChans:length(header.' ...
                  headerLeadStr{iter} '));']);
          end

          for i = 2:numChans
              if ~isequal(headSame(1, :), headSame(i,:))
                  error('Unequal top header elements among channels. Aborting.');
              end
          end
          
          %Do this for each top header element to save processing
          eval(['tempHeader(iter,:)= header.' headerLeadStr{iter} '(1, 1:numChans:length(header.' ...
              headerLeadStr{iter} '));']);
          eval(['header.' headerLeadStr{iter} '= tempHeader(iter,:);']);
      end
  end

  %VI120910A: Detect/handle header-only operation (don't read data)
  if nargout <=1
      return;
  end

  %% Process Flags

  %Determine channels to extract
  if any(ismember({'channel' 'channels'},flags))
      selectedChans = getArg({'channel' 'channels'});
      
      if ~isempty(setdiff(selectedChans,savedChans))
          selectedChans(find(setdiff(selectedChans,savedChans))) = [];
          warning('Some specified channels to extract not detected in file and, hence, ignored');
          if isempty(selectedChans)
              warning('No saved channels are specified to extract. Aborting.');
              return;
          end
      end
  else
      selectedChans = savedChans;
  end

  %This mode stays given the nature of non-selected channel storage
  %Auxiliary mapping for channel selection to index
  chanKey = num2cell(savedChans);
  chanVal = 1:length(savedChans);   %+++ Change to savedChans for selection if no resizing occurs?
  chanMap = containers.Map(chanKey,chanVal);

  %Determine slices to extract
  if numSlices >= 1 && any(ismember({'slice' 'slices'},flags))
      selectedSlices = selectImages({'slice' 'slices'},numSlices);
  else
      %Extract all slices
      selectedSlices = 1:numSlices;
  end

  % RRR Extract all frames for now
  %Determine frames to extract
  if numFrames >= 1 && any(ismember({'frame' 'frames'},flags))
      selectedFrames = selectImages({'frame' 'frames'},numFrames);
  else
      %Extract all frames
      selectedFrames = 1:numFrames;
  end


  %Determine volumes to extract
  if numVolumes >= 1 && any(ismember({'volume' 'volumes'},flags))
      selectedVolumes = selectImages({'volume' 'volumes'},numVolumes);
  else
      %Extract all frames
      selectedVolumes = 1:numVolumes;
  end

      function selection = selectImages(selectionFlags, numItems)
          if any(ismember(selectionFlags,flags))
              selection = getArg(selectionFlags);
              %Handle 'inf' specifier in slice array
              if find(isinf(selection))
                  selection(isinf(selection)) = [];
                  if max(selection) < numItems
                      selection = [selection (max(selection)+1):numItems];
                  end
              end
              if max(selection) > numItems
                  error('Frame, slice or volume values specified are not found in file');
              end
          else
              selection = 1:numItems;
          end
      end

  %Determine if any selection is being made
  forceSelection = any(ismember({'channel' 'channels' 'slice' 'slices' 'frame' 'frames' 'volume' 'volumes'},flags));

  %% Preallocate image data
  switch hTif.getTag('SampleFormat')
      case 1
          imageDataType = 'uint16';
      case 2
          imageDataType = 'int16';
      otherwise
          assert('Unrecognized or unsupported SampleFormat tag found');
  end

  %Look-up values for faster operation
  lenSelectedFrames = length(selectedFrames);
  lenSelectedChans = length(selectedChans);
  lenSelectedSlices = length(selectedSlices);
  lenSelectedVolumes = length(selectedVolumes);

  lenTotalChans = length(savedChans);
  lenTotalSlices = numSlices;
  lenTotalFrames = numFrames;
  % lenTotalVolumes = numVolumes;

  %HACK! For now there seems to be an issue with the flyback possibly due to mroi
  %still being developed. We need to take only the last section of the following values: 
  %The following also takes care of MROI mode discrepancies, since we don't have access
  %to the properties of MROI captures through the TIFF header at the moment
  numLines = hTif.getTag('ImageLength');
  numPixels = hTif.getTag('ImageWidth');

  Aout = zeros(numLines,numPixels,lenSelectedChans,lenSelectedFrames,lenSelectedSlices,lenSelectedVolumes,imageDataType);    

  %% Read image data
  selectedChans = selectedChans';

  if streamOutput
    % This mode is for the case in which the selection parameters cannot be 
    % trusted. For instance, when the number of images is different than 
    % expected, but we would still like to 
    % Checking this mode has priority given that it will always output existing data
    % No postprocessing for data (such as removing discard frames) at this point
    most.idioms.warn('Insufficient or incorrect header data.')

    %% Preallocate image data
    Aout = zeros(numLines,numPixels,numImages,imageDataType);    

    for idx = 1:numImages
      hTif.setDirectory(idx);
      Aout(:,:,idx) = hTif.read();
    end

    most.idioms.warn('Returning default, uncategorized stream of Tiff frames')

  elseif forceSelection
      for p = 1:lenSelectedVolumes
          for j = 1:lenSelectedSlices
              for k = 1:lenSelectedFrames
                  for i = 1:lenSelectedChans
                      %SELECTION MODE: (can allow parameter selection)
                      idx = chanMap(selectedChans(i));
                      %Get the tiff-index for the frames
                      idx = lenTotalChans*(selectedFrames(k) - 1) + idx;
                      %Get the tiff-index for the slices
                      idx = lenTotalFrames*lenTotalChans*(selectedSlices(j) - 1) + idx;
                      %Get the tiff-index for the volumes
                      idx = lenTotalSlices*lenTotalFrames*lenTotalChans*(selectedVolumes(p) - 1) + idx;
                      
                      %+++ Test the following expression.
                      if ismember(selectedChans(i), savedChans)
                          hTif.setDirectory(idx);
                          Aout(:,:,i,k,j,p) = hTif.read();
                      end
                  end
              end
          end
      end
  else
      idx = 0;
      for p = 1:lenSelectedVolumes
          for j = 1:lenSelectedSlices
              for k = 1:lenSelectedFrames
                  for i = 1:lenSelectedChans
                      %NO-SELECTION MODE: (more efficient)
                      idx = idx + 1;

                      if ismember(selectedChans(i), savedChans)
                          hTif.setDirectory(idx);
                          Aout(:,:,i,k,j,p) = hTif.read();
                      end
                  end
              end
          end
      end
  end

  % Prepare imgInfo
  imgInfo.numImages = numImages;
  imgInfo.numChans = numChans;
  imgInfo.numPixels = numPixels;
  imgInfo.numLines = numLines;
  imgInfo.numSlices = numSlices;
  imgInfo.numVolumes = numVolumes;
  imgInfo.numFrames = numFrames;
  imgInfo.filename = fileName;	
  imgInfo.si_ver = si_ver;	
        
  %% GENERAL HELPERS
  function arg = getArg(flag)
    [tf,loc] = ismember(flag,flags); %Use this approach, instead of intersect, to allow detection of flag duplication
    if length(find(tf)) > 1
      error(['Flag ''' flag ''' appears more than once, which is not allowed']);
    else %Extract location of specified flag amongst flags
      loc(~loc) = [];
    end
    flagIndex = flagIndices(loc);
    if length(varargin) <= flagIndex
      arg = [];
      return;
    else
      arg = varargin{flagIndex+1};
      if ischar(arg) && ismember(lower(arg),flags) %Handle case where argument was omitted, and next argument is a flag
        arg = [];
      end
    end
  end

  function [s, si_ver, numhdrs, headerLeadStr, numImg] = parseHeaderToStruct(tifObj)
    s = struct();
    numhdrs = 0;
    headerLeadStr = {};
    numImg = 0;
    si_ver = 0;
		auxTriggerCell = cell(0);

    % Before anything else, see if the tiff file has any image-data
    try
      %Parse SI from the first frame
      numImg = 1;
      while ~tifObj.lastDirectory()
        tifObj.nextDirectory();
        numImg = numImg + 1;
      end
      tifObj.setDirectory(1);
    catch
      warning('The tiff file may be corrupt.')
    end

    %Make sure the tiff file's ImageDescription didn't go over the limit set in 
    %Acquisition.m:LOG_TIFF_HEADER_EXPANSION
    try
      if ~isempty(strfind(tifObj.getTag('ImageDescription'), '<output truncated>'))
        most.idioms.warn('Corrupt header data');
        return
      end
    catch
      most.idioms.warn('Corrupt or incomplete tiff header');
      return
    end

    % Then check if the header is empty data
    try
      frameString= tifObj.getTag('ImageDescription');
    catch 
       most.idioms.warn('The input tiff file may be corrupt or its header empty')
       return;
    end
    
    si_ver = getHeaderProperty(frameString,'scanimage.SI.VERSION_MAJOR');
    if isempty(si_ver)
      si_ver = getHeaderProperty(frameString,'scanimage.SI5.VERSION_MAJOR');
    end

    if isempty(si_ver)
      most.idioms.warn('This file version is not supported by opentif. Returning available information.');
      si_ver = 0;
      return;
    elseif si_ver == 5
      headerLeadStr = {   'frameNumbers',...
                          'frameTimestamps',...
                          'acqTriggerTimestamps',...
                          'nextFileMarkerTimestamps',...
                          'dcOverVoltage'};
      numhdrs = length(headerLeadStr);
    elseif si_ver == 2015
      headerLeadStr = {   'frameNumbers',...
                          'acquisitionNumbers',...
                          'frameNumberAcquisition',...
                          'frameTimestamps',...
                          'acqTriggerTimestamps',...
                          'nextFileMarkerTimestamps',...
                          'endOfAcquisition',... 
                          'endOfAcquisitionMode',...
                          'dcOverVoltage'};

      si2015_scanMode = getHeaderProperty(frameString,'scanimage.SI.imagingSystem');
      if strcmp(si2015_scanMode, 'Resonant')
        auxTriggerCell = { 	'auxTrigger0',...
                  'auxTrigger1',...
                  'auxTrigger2',...
                  'auxTrigger3'};
      end
    end

    numhdrs = length(headerLeadStr);

    rows = textscan(frameString,'%s','Delimiter','\n');            
    rows = rows{1};

    %If the first frame is empty return
    if isempty(rows)
      return;
    end
    
    for c = 1:numhdrs
      % Unassigned values are set to zero
      eval(['s.' headerLeadStr{c} '=zeros(1,numImg);'])
    end

    try
      for frame = 1:numImg
        frameString  = tifObj.getTag('ImageDescription');       
        rows = textscan(frameString,'%s','Delimiter','\n');            
        rows = rows{1};

        for c = 1:numhdrs
          row = rows{c};
          % replace top-level name with 'obj'
          [~, rmn] = strtok(row,'=');

          % Check if there is a value to assign
          if strcmp(strtrim(rmn),'=') ;
            % This unassigned parameter value will be set to 0
            continue;
          else
            evalStr = ['s.' headerLeadStr{c} '(' num2str(frame) ')'  rmn];
          end
          eval([evalStr ';']);
        end
        
        % Process auxTriggers if any
        if si_ver == 2015 && strcmp(si2015_scanMode, 'Resonant') && numel(auxTriggerCell) == 4
            for c = 1:numel(auxTriggerCell)
                row = rows{c+numhdrs};

                % replace top-level name with 'obj'
                [~, rmn] = strtok(row,'=');
                
                % Check if there is a value to assign
                if strcmp(strtrim(rmn),'=') ;
                    % This unassigned parameter value will be set to 0
                    continue;
                else
                    evalStr = ['s.' auxTriggerCell{c} '{' num2str(frame) '}'  rmn];
                end
                eval([evalStr ';']);
            end
        end

        if frame ~= numImg
          tifObj.nextDirectory();
        end
      end
    catch err
      % Give more information for mismatch.
      most.idioms.warn('The file has an unexpected header format');
      %most.idioms.warn(err.message);
      %Clear struct since format is unexpected
      %+++ Find a better way to handle this
      s = struct();
      return
    end  % end try/catch
          
    % Handle SI field            
    tifObj.setDirectory(1);
    frameString  = tifObj.getTag('ImageDescription');       
    rows = textscan(frameString,'%s','Delimiter','\n');            
    rows = rows{1};
    
    for c = numhdrs+1:numel(rows)
      row = rows{c};
      
      % replace top-level name with 'obj'
      [~, rmn] = strtok(row,'.');
      row = ['s' rmn];

      % deal with nonscalar nested structs/objs
      pat = '([\w]+)__([0123456789]+)\.';
      replc = '$1($2).';
      row = regexprep(row,pat,replc);

      % handle unencodeable value or nonscalar struct/obj
      unencodeval = '<unencodeable value>';
      if strfind(row,unencodeval)
        row = strrep(row,unencodeval,'[]');
      end
      nonscalarstructobjstr = '<nonscalar struct/object>';
      if strfind(row,nonscalarstructobjstr)
        row = strrep(row,nonscalarstructobjstr,'[]');
      end

      % handle ND array format produced by array2Str
      try 
        if ~isempty(strfind(row,'&'))
          equalsIdx = strfind(row,'=');
          [dimArr rmn] = strtok(row(equalsIdx+1:end),'&');
          arr = strtok(rmn,'&');
          arr = reshape(str2num(arr),str2num(dimArr)); %#ok<NASGU,ST2NM>
          eval([row(1:equalsIdx+1) 'arr;']);
        else
          eval([row ';']);
        end
      catch ME %Warn if assignments to no-longer-extant properties are found
        if strcmpi(ME.identifier,'MATLAB:noPublicFieldForClass')
          equalsIdx = strfind(row,'=');
          warnMsg = sprintf(1,'Property ''%s'' was specified, but does not exist for class ''%s''\n', deblank(row(3:equalsIdx-1)),class(s));
          most.idioms.warn(warnMsg);
        else
          ME.rethrow();
        end
      end
    end
  end

  function s = extractHeaderData(header, si_version)
    if isfield(header,'SI')
      localHdr = header.SI;
    elseif isfield(header,'SI5')
      localHdr = header.SI5;
    else
      assert(false);
    end

    if si_version == 5
      s.savedChans = localHdr.channelsSave;
      s.numPixels = localHdr.pixelsPerLine;
      s.numLines = localHdr.linesPerFrame;
      s.numVolumes = localHdr.fastZNumVolumes;
      if isfield(localHdr,'acqNumAveragedFrames')
        saveAverageFactor = localHdr.acqNumAveragedFrames;
      else
        assert(false);
      end
      s.numFrames = localHdr.acqNumFrames / saveAverageFactor;
      s.numSlices = localHdr.stackNumSlices;

    elseif si_version == 2015
      %+++ FIX ME
      s.savedChans = localHdr.hChannels.channelSave;
      s.numPixels = localHdr.hRoiManager.pixelsPerLine;
      s.numLines = localHdr.hRoiManager.linesPerFrame;

      if localHdr.hFastZ.enable
        s.numVolumes = localHdr.hFastZ.numVolumes;
        s.numSlices = localHdr.hStackManager.slicesPerAcq;
        s.numFrames = 1;

        % Assuming that we only have discard frames during FastZ acquisitions
        s.discardFlybackframesEnabled = localHdr.hFastZ.discardFlybackFrames;
        s.numDiscardFrames = localHdr.hFastZ.numDiscardFlybackFrames; 
        s.numFramesPerVolume = localHdr.hFastZ.numFramesPerVolume;  %Includes flyback frames
      else
        s.numVolumes = 1;
        s.numFrames = localHdr.hStackManager.framesPerSlice;
        s.numSlices = localHdr.hStackManager.slicesPerAcq;

        s.discardFlybackframesEnabled = false;
        s.numDiscardFrames = localHdr.hFastZ.numDiscardFlybackFrames;    
        s.numFramesPerVolume = localHdr.hFastZ.numFramesPerVolume;  %Includes flyback frames
      end

      % NOTE: This assumes you are using tiff files generated on non-simulated
      %       mode. In this case, non-FastZ tiff files seem to differ between these modes
      if s.numSlices > 1
          s.numFrames = 1;
      end


    else
      assert(false);
    end
  end

  function forceQuit = streamOutputQuit()
    %This function returns available data and promts the program to exit
    %+++ output parameter should be generalized
    forceQuit = true;

    %% Preallocate image data
    switch hTif.getTag('SampleFormat')
      case 1
        imageDataType = 'uint16';
      case 2
        imageDataType = 'int16';
      otherwise
        assert('Unrecognized or unsupported SampleFormat tag found');
    end

    numLines = hTif.getTag('ImageLength');
    numPixels = hTif.getTag('ImageWidth');

    Aout = zeros(numLines,numPixels,numImages,imageDataType);    
    imgInfo.numImages = numImages;	% Only the number of images is reliable
    imgInfo.filename = fileName;	% As well as the filename, of course
    imgInfo.si_ver = si_ver;	% ScanImage version 

    for idx = 1:numImages
      hTif.setDirectory(idx);
      Aout(:,:,idx) = hTif.read();
    end



    most.idioms.warn('Returning default, uncategorized stream of Tiff frames')
    return
  end

  function val = getHeaderProperty(imdescription,propfullname)
    try
      str = regexpi(imdescription,sprintf('(?<=%s ?= ?).*$',propfullname),'match','once','lineanchors','dotexceptnewline');
    catch
      str = [''''';'];
    end
    if isempty(str);
      str = [''''';'];
    end
    val = eval(str);
  end

end


%--------------------------------------------------------------------------%
% opentif.m                                                                %
% Copyright © 2015 Vidrio Technologies, LLC                                %
%                                                                          %
% ScanImage 2015 is premium software to be used under the purchased terms  %
% Code may be modified, but not redistributed without the permission       %
% of Vidrio Technologies, LLC                                              %
%--------------------------------------------------------------------------%

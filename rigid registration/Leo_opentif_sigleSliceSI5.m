function [hTif, header,Aout] = Leo_opentif_sigleSliceSI5(varargin)
%% function [header,Aout] = open_tif(varargin)
% Opens a ScanImage TIF file, extracting its header information and, if specified, stores all of image contents as output array Aout if specified. 
% By default, Aout, if specified for output, is of size MxNxCxK,where C spans the channel indices, and K the slice/frame indices.
%
% NOTE: IF the second output argument (Aout) is not assigned to output variable
%       THEN image file is not actually read -- only  header information is extracted
%
%% SYNTAX
%   open_tif()
%   open_tif(filename)
%   header = open_tif(...)
%   [header,Aout] = open_tif(...)
%       filename: Name of TIF file, with or without '.tif' extension. If omitted, a dialog is launched to allow interactive selection.
%       flagN/flagNArg: Flags (string-valued) and/or flag/value pairs, in any order, specifying options to use in opening specified file
%
%       header: Structure comprising information stored by ScanImage into TIF header
%       Aout: MxNxCxK array, with images of size MxN for each of C colors and K slices or frames. Default type is uint16.
%
% NOTE: IF the second output argument (Aout) is not assigned to output variable
%       THEN image file is not actually read -- only header information is extracted
%
%% FLAGS (case-insensitive)
%
%   WITH ARGUMENTS
%       'channel' or 'channels': Argument specifies subset of channel(s) to extract. Ex: 1,[1 3], 2:4. 
%       'slice' or 'slices': Argument specifies subset of slices present to extract. Use 'inf' to specify all slices above highest specified value. Ex: 1:30, [50 inf], [1:9 11:19 21 inf]
%       'frame' or 'frames': Argument specifies subset of frames present to extract. Use 'inf' to specify all frames above highest specified value. Ex: 1:30, [50 inf], [1:9 11:19 21 inf]
%
%% NOTES
%   This function replaces the scim_openTif() function supplied with ScanImage 4.2
%
%   TODO: Port more advaced features to ScanImage 5
%   TODO: Add support for ScanImage 3 and 4 tif files
%   TODO: Add output style flags for Aout
%   TODO: Add output style flags
%
%% CREDITS
%   Based on scim_openTif on 4/16/09, by Vijay Iyer
%
%% Constants/Inits
error(nargoutchk(0,3,nargout,'struct'));

%% Parse input arguments
flagNames = {'channel' 'channels' 'slice' 'slices' 'frame' 'frames'};
argFlags = {'channel' 'channels' 'slice' 'slices' 'frame' 'frames'};

flagIndices = find(cellfun(@(x)ischar(x) && (ismember(lower(x),flagNames) || ismember(lower(x),argFlags)),varargin));
flags = cellfun(@lower,varargin(flagIndices),'UniformOutput',false);
if isempty(flags)
    flags = {};
end

%% Determine input file
if isempty(find(flagIndices==1)) && nargin>=1 && ischar(varargin{1})
    fileName = varargin{1};
else
    fileName = '';
end

if isempty(fileName)
    [f, p] = uigetfile({'*.tif;*.tiff'},'Select Image File');
    if f == 0
        return;
    end
    fileName = fullfile(p,f); 
end

%Extract filepath for future use
[filePath,fileStem,fileExt] = fileparts((fileName));

%% Read TIFF file; extract # frames & image header
if ~exist(fileName,'file')
    error('''%s'' is not a recognized flag or filename. Aborting.',fileName);
end

warning('off','MATLAB:tifflib:TIFFReadDirectory:libraryWarning');
hTif = Tiff(fileName);

fileVersion = 5;
header = parseHeaderToStruct(hTif,fileVersion);            
hdr = extractHeaderData(header,fileVersion);
 
if nargout <=1
    return;
end
%% Read image meta-data
savedChans = hdr.savedChans;
numChans = length(savedChans);
numPixels = hdr.numPixels;
numLines = hdr.numLines;
numSlices = hdr.numSlices;
numFramesPerSlice = hdr.numFrames;
%numTotalFrames = numSlices*numFramesPerSlice;
numTotalFrames = numSlices*length( imfinfo(fileName));%LTP MODIFIED

if ~numFramesPerSlice || ~numSlices
    error('Acquisition did not complete a single frame or slice. Aborting.');
end

%% Process Flags

if any(ismember({'slice' 'slices'},flags)) && any(ismember({'frame' 'frames'},flags))
    error('Simultaneous queries of slices and frames not supported at this time')
end

%Determine channels to extract
if any(ismember({'channel' 'channels'},flags))
    selChans = getArg({'channel' 'channels'});
    
    if ~isempty(setdiff(selChans,savedChans))
        selChans(find(setdiff(selChans,savedChans))) = [];
        warning('Some specified channels to extract not detected in file and, hence, ignored');
        if isempty(selChans)
            warning('No saved channels are specified to extract. Aborting.');
            return;
        end
    end
    reducedChans = length(selChans) < length(savedChans) %Determine if # channels was reduced by  
else
    selChans = savedChans;
end

%Determine slices & frames to extract
if numSlices >= 1 && any(ismember({'slice' 'slices'},flags))
    preSelection = selectImages({'slice' 'slices'},numSlices);
    selection = [];
    for i = 1:length(preSelection)
        selection = [selection (preSelection(i)-1)*numFramesPerSlice+1:...
            (preSelection(i))*numFramesPerSlice];
    end
    selection = sort(selection);
%     selectionStr = 'Slice';
elseif numTotalFrames >= 1 && any(ismember({'frame' 'frames'},flags))
    selection = selectImages({'frame' 'frames'},numTotalFrames);
%     selectionStr = 'Frame';
else
    %Extract all frames
    selection = 1:numTotalFrames;
%     selectionStr = '';
end

numSelections = length(selection);
    
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
                error('Slice or frame values specified are not found in file');
            end
        else
            selection = 1:numItems;
        end
    end

%% Preallocate image data
switch hTif.getTag('SampleFormat')
    case 1
        imageDataType = 'uint16';
    case 2
        imageDataType = 'int16';
    otherwise
        assert('Unrecognized or unsupported SampleFormat tag found');
end

Aout = zeros(numLines,numPixels,length(selChans),numSelections,imageDataType);

%% Read image data
for i = 1:length(selection)
    for j = 1:length(selChans)
        idx = numChans * (selection(i) - 1) + selChans(j);
        if ismember(selChans(j), savedChans) %change?
            hTif.setDirectory(idx);
            Aout(:,:,j,i) = hTif.read();
        end
    end
end

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

    function s=parseHeaderToStruct(tifObj,fileVersion)
        if fileVersion == 5
            s = struct();
    
            %Parse SI5 from the first frame
            numImages = 1;
            while ~tifObj.lastDirectory()
                numImages = numImages + 1;
                tifObj.nextDirectory();
            end
            tifObj.setDirectory(1);
            
            frameString  = tifObj.getTag('ImageDescription');       
            rows = textscan(frameString,'%s','Delimiter','\n');            
            rows = rows{1};

            %If the first frame is empty return
            if isempty(rows)
                return;
            end
            
            %Parse the 5 array vars corresponding to the first 5 entries
            headerLead = {  'frameNumbers',...
                            'frameTimestamps',...
                            'acqTriggerTimestamps',...
                            'nextFileMarkerTimestamps',...
                            'dcOverVoltage'};
                            
            %Loop for all frames
            for c = 1:5
                eval(['s.' headerLead{c} '=zeros(1,numImages);'])
            end
            for frame = 1:numImages
                frameString  = tifObj.getTag('ImageDescription');       
                rows = textscan(frameString,'%s','Delimiter','\n');            
                rows = rows{1};
                
                for c = 1:5
                    row = rows{c};

                    % replace top-level name with 'obj'
                    [~, rmn] = strtok(row,'=');
                    row = ['s.' headerLead{c} '(frame)'  rmn];
                    %Check for empty cases +++
                    eval([row ';']);
                end
                if frame ~= numImages
                    tifObj.nextDirectory();
                end
            end
            
            % Handle SI5 field            
            
            tifObj.setDirectory(1);
            frameString  = tifObj.getTag('ImageDescription');       
            rows = textscan(frameString,'%s','Delimiter','\n');            
            rows = rows{1};
            
            for c = 6:numel(rows)
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
                        fprintf(1,'WARNING: Property ''%s'' was specified, but does not exist for class ''%s''\n', deblank(row(3:equalsIdx-1)),class(s));
                    else
                        ME.rethrow();
                    end
                end
            end
        end
    end

    function s = extractHeaderData(header,fileVersion)
        
        if fileVersion == 5
            if isfield(header,'SI5')
                localHdr = header.SI5;
            else
                assert(false);
            end
            
            s.savedChans = localHdr.channelsSave;
            s.numPixels = localHdr.pixelsPerLine;
            s.numLines = localHdr.linesPerFrame;
            
            if isfield(localHdr,'acqNumAveragedFrames')
                saveAverageFactor = localHdr.acqNumAveragedFrames;
            else
                assert(false);
            end

            s.numFrames = localHdr.acqNumFrames / saveAverageFactor;
            
            s.numSlices = localHdr.stackNumSlices;
            
        else
            assert(false);
        end

    end


end
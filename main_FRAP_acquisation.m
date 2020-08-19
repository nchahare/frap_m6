% File directory for analysing FLIP data from mic6 with volume correction 
% equations. Check the math.pdf for the math.
% Created on 
% 13-Aug-2020 17:19:55
% @author: nchahare
%  
% main
% ---display image
% ---mark contours around the image
% ---extracting data
% ---calucating parameters from the equations
% ---storing the data

clear all; clc; close all

%-------------------- START OF THE CODE -----------------------------------

%% <<<<<<<<<<<<<<<<<<< set the directory

% record current location and then select the directory. change directory
% to that one and record all the filenames and location. Then come back to
% current location
currentdirectory = pwd;
% directory =...
%     'C:\Users\NChahare\OneDrive - IBEC\Lab Stuff\FLIP\FLIP2020\equations adjusted for volume\data\frap flip example 2020_07_21\FRAP'
directory = uigetdir('MATLAB Root Folder');
cd(directory);
Files = dir('*.czi');
cd(currentdirectory)


% domore means analyse more and imagenumber refers to individual file. I
% have limited to do it till 10 files but you can increase it to 50 if you
% want
for imagenumber = 1: length(Files)

    % asking you to continue or quit
    promptMessage = sprintf('Do you want to analyse file number %d,\nor Quit?', imagenumber);
    button = questdlg(promptMessage, 'Continue?', 'Quit', 'Continue','Quit');
    if strcmpi(button, 'Quit')
        break;
    end
    



    %% <<<<<<<<<<<<<<<<<<< reading file
    % reading the individual file
    File.CZI = fullfile(Files(imagenumber).folder,Files(imagenumber).name);
    File.TIF = File.CZI; File.TIF(end-2:end) = 'tif';

    % reading xml data
    InfoCZI  = czifinfo(File.CZI);
    InfoCZI  = xml2struct(InfoCZI.metadataXML);

    %% <<<<<<<<<<<<<<<<<<< openning tif stack
    InfoImage    = imfinfo(File.TIF);    %getting image info
    NumberImages = length(InfoImage);    %getting total count of images

    % saving all the images in 'Images'
    Images  = zeros(InfoImage(1).Height,InfoImage(1).Width,...
    NumberImages,'uint16');
    TifLink = Tiff(File.TIF, 'r');

    for i=1:NumberImages
        TifLink.setDirectory(i);
        Images(:,:,i)=TifLink.read();
    end
    TifLink.close();
    clear i TifLink

    % finding imaging interval
    Imagejinfo = InfoImage(1).ImageDescription;
    kk = strfind(Imagejinfo,'finterval=');


    %% <<<<<<<<<<<<<<<<<<< reading data about layers
    % A lot of work had to be done for this to find each layer and type. i am
    % sick of this, glad that it is done. DONOT MESS WITH THIS!!! Skip this
    % section
    Layers = InfoCZI.ImageDocument.Metadata.Layers.Layer;

    flag = 0;
    for i = 1:length(Layers)
        A =  fieldnames(Layers{1,i}.Elements);
        for j = 1:length(A)
            Z = getfield(Layers{1,i}.Elements,A{j});
            for k = 1:length(Z)
                flag = flag +1;
                Elements(flag).ROItype = A{j};
                if iscell(Z)
                    Elements(flag).Property = getfield(Z{k},'Geometry');
                else
                    Elements(flag).Property = getfield(Z(k),'Geometry');
                end
            end
        end
    end
    clear Z i j k flag Layers A

    % finding the frap circle radius and center
    for i = 1:length(Elements)
        Radius(i) = str2num(Elements(i).Property.Radius.Text);
    end
    [frap.Radius indofElements] = max(Radius);
    frap.Center = [str2num(Elements(indofElements).Property.CenterX.Text),...
                   str2num(Elements(indofElements).Property.CenterY.Text)];


    %% <<<<<<<<<<<<<<<<<<< mark cell boundary as an ROI
    % showing the image in jet color and showing original frap zone and
    % interest zone with circles
    figure
    [rows, columns, ~] = size(Images(:,:,1));
    imagesc(Images(:,:,1)); axis equal; axis off; colormap(jet);
    drawcircle('Center',frap.Center,...
    'Radius', frap.Radius,'Color','r'); % drawing the circle
    drawcircle('Center',frap.Center,...
    'Radius', 3*frap.Radius,'Color','b'); % drawing the circle
    title({'Original Image. Draw the cell boundary';...
    'Then double click in the middle to accept it.'},...
    'FontSize', 18);
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % Maximize figure.
    set(gcf,'name','Image Analysis Demo','numbertitle','off')

    % Ask user to draw polygons
    % Create a binary image for all the regions we will draw.
    cellmask = false(rows, columns);
    [thisSinglePolygonImage, xi, yi] = roipoly();
    cellmask = cellmask | thisSinglePolygonImage;
    cellmask = uint16(cellmask);
    close all

    % % showing mask for the cell
    % figure; imshow(cellmask);

    clear xi yi thisSinglePolygonImage originalImage...
    caption button again titleBarCaption hc k


    %% <<<<<<<<<<<<<<<<<<< Frap region mask
    figure; % applying mask for the first image
    imagesc(Images(:,:,1).*cellmask); axis equal; axis off; colormap(jet);
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % Maximize figure.
    title({'check the ROI'})
    hc = drawcircle('Center',frap.Center,...
    'Radius', frap.Radius,'Color','r'); % drawing the circle
    frapMask = uint16(createMask(hc));       % create a mask for frapzone
    % figure; imagesc(frapMask); axis equal; axis off;

    clear hc 


    %% <<<<<<<<<<<<<<<<<<< storing relevant data in RAWDATA structure
    RAWDATA(imagenumber).File           = File.CZI;     % file name and location
    RAWDATA(imagenumber).Images         = Images;       % Raw images in a 3D array
    RAWDATA(imagenumber).CellMask       = cellmask;     % cell mask
    RAWDATA(imagenumber).FrapMask       = frapMask;     % frap mask
    RAWDATA(imagenumber).FrapRadius     = frap.Radius;  % radius of frap
    RAWDATA(imagenumber).FrapCenter     = frap.Center;  % centre of frap
    RAWDATA(imagenumber).xscale         = ...           % scale of image
    str2num(InfoCZI.ImageDocument.Metadata.Experiment.ExperimentBlocks. ...
    AcquisitionBlock.AcquisitionModeSetup.ScalingX.Text);
    RAWDATA(imagenumber).yscale = ...
    str2num(InfoCZI.ImageDocument.Metadata.Experiment.ExperimentBlocks. ...
    AcquisitionBlock.AcquisitionModeSetup.ScalingY.Text);
    RAWDATA(imagenumber).timeinterval =...              % time interval
    str2num(Imagejinfo(kk+10:kk+13));

    clear ans cellmaks Elements File frap cellmask i Imagejinfo Images ...
    indofElements InfoCZI InfoImage kk mask nucleolimask NumberImages ...
    promptMessage Radius


    %% <<<<<<<<<<<<><<<<<< collecting values for intensity
    FinalImage = RAWDATA(imagenumber).Images;    
    cellmask = RAWDATA(imagenumber).CellMask;
    frapmask = RAWDATA(imagenumber).FrapMask;
    areascale = RAWDATA(imagenumber).xscale^2;

    for i = 1: size(FinalImage,3)
        MaskedImage = FinalImage(:,:,i).*(cellmask); % applying mask
        Fraped = MaskedImage.*frapmask;              % applying frap mask
        Frapedarray = Fraped(logical(frapmask));     % all the info in the frap mask
        Cellarray = MaskedImage(logical(cellmask));  % all the info in the cell mask

        FrapIntensity(i,1).time = i*RAWDATA(imagenumber).timeinterval;
        FrapIntensity(i,1).mean  = mean(Frapedarray);
        FrapIntensity(i,1).median  = median(Frapedarray);
        FrapIntensity(i,1).IntDen  = mean(Frapedarray)*length(Frapedarray)*areascale;

        CellIntensity(i,1).time = FrapIntensity(i,1).time;
        CellIntensity(i,1).mean  = mean(Cellarray);
        CellIntensity(i,1).median  = median(Cellarray);
        CellIntensity(i,1).IntDen  = mean(Cellarray)*length(Cellarray)*areascale;
    end

    RAWDATA(imagenumber).FrapIntensity = FrapIntensity;
    RAWDATA(imagenumber).CellIntensity = CellIntensity;

    clear areascale  Cellarray CellIntensity FinalImage Fraped Frapedarray...
    FrapIntensity frapmask frapMask MaskedImage i


    %% <<<<<<<<<<<<><<<<<< data quality check
    promptMessage = sprintf('Did you think that this data was good?');
    button = questdlg(promptMessage, 'Really?', 'Good', 'Shit','Good');
    RAWDATA(imagenumber).quality = button;


end

close all

% save file as user input name
prompt = {'Save data file as:'};
dlgtitle = 'Input';
dims = [1 35];
filename = [char(inputdlg(prompt,dlgtitle,dims)) '.mat'];
save(filename,'RAWDATA')


%-------------------- END OF THE CODE -------------------------------------------------------------------


%% --------------------- FUNCTIONS ----------------------------------------
% I STRONGLY RECOMMEND TO NOT READ THIS STUFF. IT IS NOT COMMENTED WELL
%% xml to structure converter

function  outStruct  = xml2struct(input)
%XML2STRUCT converts xml file into a MATLAB structure
%
% outStruct = xml2struct2(input)
% 
% xml2struct2 takes either a java xml object, an xml file, or a string in
% xml format as input and returns a parsed xml tree in structure. 
% 
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
%
% Originally written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increase by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
% Modified by X. Mo, University of Wisconsin, 12-5-2012
% Modified by Chao-Yuan Yeh, August 2016

errorMsg = ['%s is not in a supported format.\n\nInput has to be',...
        ' a java xml object, an xml file, or a string in xml format.'];

% check if input is a java xml object
if isa(input, 'org.apache.xerces.dom.DeferredDocumentImpl') ||...
        isa(input, 'org.apache.xerces.dom.DeferredElementImpl')
    xDoc = input;
else
    try 
        if exist(input, 'file') == 2
            xDoc = xmlread(input);
        else
            try
                xDoc = xmlFromString(input);
            catch
                error(errorMsg, inputname(1));
            end
        end
    catch ME
        if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
            error(errorMsg, inputname(1));
        else
            rethrow(ME)
        end
    end
end

% parse xDoc into a MATLAB structure
outStruct = parseChildNodes(xDoc);
    
end

% ----- Local function parseChildNodes -----
function [children, ptext, textflag] = parseChildNodes(theNode)
% Recurse over node children.
children = struct;
ptext = struct; 
textflag = 'Text';

if hasChildNodes(theNode)
    childNodes = getChildNodes(theNode);
    numChildNodes = getLength(childNodes);

    for count = 1:numChildNodes

        theChild = item(childNodes,count-1);
        [text, name, attr, childs, textflag] = getNodeData(theChild);
        
        if ~strcmp(name,'#text') && ~strcmp(name,'#comment') && ...
                ~strcmp(name,'#cdata_dash_section')
            % XML allows the same elements to be defined multiple times,
            % put each in a different cell
            if (isfield(children,name))
                if (~iscell(children.(name)))
                    % put existsing element into cell format
                    children.(name) = {children.(name)};
                end
                index = length(children.(name))+1;
                % add new element
                children.(name){index} = childs;
                
                textfields = fieldnames(text);
                if ~isempty(textfields)
                    for ii = 1:length(textfields)
                        children.(name){index}.(textfields{ii}) = ...
                            text.(textfields{ii});
                    end
                end
                if(~isempty(attr)) 
                    children.(name){index}.('Attributes') = attr; 
                end
            else
                % add previously unknown (new) element to the structure
                
                children.(name) = childs;
                
                % add text data ( ptext returned by child node )
                textfields = fieldnames(text);
                if ~isempty(textfields)
                    for ii = 1:length(textfields)
                        children.(name).(textfields{ii}) = text.(textfields{ii});
                    end
                end

                if(~isempty(attr)) 
                    children.(name).('Attributes') = attr; 
                end
            end
        else
            ptextflag = 'Text';
            if (strcmp(name, '#cdata_dash_section'))
                ptextflag = 'CDATA';
            elseif (strcmp(name, '#comment'))
                ptextflag = 'Comment';
            end

            % this is the text in an element (i.e., the parentNode) 
            if (~isempty(regexprep(text.(textflag),'[\s]*','')))
                if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
                    ptext.(ptextflag) = text.(textflag);
                else
                    % This is what happens when document is like this:
                    % <element>Text <!--Comment--> More text</element>
                    %
                    % text will be appended to existing ptext
                    ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
                end
            end
        end

    end
end
end

% ----- Local function getNodeData -----
function [text,name,attr,childs,textflag] = getNodeData(theNode)
% Create structure of node info.

%make sure name is allowed as structure name
name = char(getNodeName(theNode));
name = strrep(name, '-', '_dash_');
name = strrep(name, ':', '_colon_');
name = strrep(name, '.', '_dot_');
name = strrep(name, '_', 'u_');

attr = parseAttributes(theNode);
if (isempty(fieldnames(attr))) 
    attr = []; 
end

%parse child nodes
[childs, text, textflag] = parseChildNodes(theNode);

% Get data from any childless nodes. This version is faster than below.
if isempty(fieldnames(childs)) && isempty(fieldnames(text))
    text.(textflag) = char(getTextContent(theNode));
end

% This alterative to the above 'if' block will also work but very slowly.
% if any(strcmp(methods(theNode),'getData'))
%   text.(textflag) = char(getData(theNode));
% end
    
end

% ----- Local function parseAttributes -----
function attributes = parseAttributes(theNode)
% Create attributes structure.
attributes = struct;
if hasAttributes(theNode)
   theAttributes = getAttributes(theNode);
   numAttributes = getLength(theAttributes);

   for count = 1:numAttributes
        % Suggestion of Adrian Wanner
        str = char(toString(item(theAttributes,count-1)));
        k = strfind(str,'='); 
        attr_name = str(1:(k(1)-1));
        attr_name = strrep(attr_name, '-', '_dash_');
        attr_name = strrep(attr_name, ':', '_colon_');
        attr_name = strrep(attr_name, '.', '_dot_');
        attributes.(attr_name) = str((k(1)+2):(end-1));
   end
end
end

% ----- Local function xmlFromString -----
function xmlroot = xmlFromString(iString)
import org.xml.sax.InputSource
import java.io.*

iSource = InputSource();
iSource.setCharacterStream(StringReader(iString));
xmlroot = xmlread(iSource);
end

%% CZI metadata reader

function  fileInfo  = czifinfo( filename, varargin )
%CZIFINFO returns informaion of Zeiss CZI file
%
%   czifinfo returns information of czi file includingl pixel type,
%   compression method, fileGUID, file version number, a structure
%   recording various information of raw image data including data start
%   position within the czi file, data size and spatial coordinates. Also
%   the function returns associated metadata in the field named 'XML_text'.
%   This can be saved as an .xml file and examined in web browser. 
%
%   Version 1.0
%   Copyright Chao-Yuan Yeh, 2016

fID = fopen(filename);

while true
    segHeader = readSegHeader(fID);
    if strfind(segHeader.SID, 'ZISRAWSUBBLOCK')
        fileInfo.genInfo = readMinSUBBLOCKHeader(fID);
        break
    end
    fseek(fID, segHeader.currPos + segHeader.allocSize, 'bof');
end
count = 0;
frewind(fID);
flag = 1;
sBlockCount_P0 = 0;
sBlockCount_P2 = 0;
while flag
    segHeader = readSegHeader(fID);
    if segHeader.allocSize
        if strfind(segHeader.SID, 'ZISRAWSUBBLOCK')
            [sBlockHeader, pyramidType] = readPartSUBBLOCKHeader(fID);
            switch pyramidType
                case 0
                    sBlockCount_P0 = sBlockCount_P0 + 1;
                    fileInfo.sBlockList_P0(sBlockCount_P0) = sBlockHeader;
                case 2
                    if strcmpi(varargin,'P2')
                        sBlockCount_P2 = sBlockCount_P2 + 1;
                        fileInfo.sBlockList_P2(sBlockCount_P2) = sBlockHeader;
                    end
            end
        elseif strfind(segHeader.SID, 'ZISRAWFILE')
            fileInfo.fileHeader = readFILEHeader(fID);
        elseif strfind(segHeader.SID, 'ZISRAWATTACH')
            count = count + 1;
            readAttach(fID);
        end
        flag = fseek(fID, segHeader.currPos + segHeader.allocSize, 'bof') + 1;
    else
        flag = 0;
    end
end

fseek(fID, 92, 'bof');
fseek(fID, fileInfo.fileHeader.mDataPos, 'bof');
fseek(fID, fileInfo.fileHeader.mDataPos + 32, 'bof');
XmlSize = uint32(fread(fID, 1, '*uint32'));
fseek(fID, fileInfo.fileHeader.mDataPos + 288, 'bof');
fileInfo.metadataXML = fread(fID, XmlSize, '*char')';
fclose(fID);
disp(count)
end

function segHeader = readSegHeader(fID)
segHeader.SID = fread(fID, 16, '*char')';
segHeader.allocSize = fread(fID, 1, '*uint64');
fseek(fID, 8, 'cof'); 
segHeader.currPos = ftell(fID);
end

function sBlockHeader = readMinSUBBLOCKHeader(fID)
fseek(fID, 18, 'cof');
sBlockHeader.pixelType = getPixType(fread(fID, 1, '*uint32'));
fseek(fID, 12, 'cof');
sBlockHeader.compression = getCompType(fread(fID, 1, '*uint32'));
fseek(fID, 6, 'cof');
sBlockHeader.dimensionCount = fread(fID, 1, '*uint32');
end

function [sBlockHeader, pyramidType] = readPartSUBBLOCKHeader(fID)
currPos = ftell(fID);
mDataSize = fread(fID, 1, '*uint32');
fseek(fID, 4, 'cof');
sBlockHeader.dataSize = fread(fID, 1, '*uint64');
fseek(fID, 22, 'cof');
pyramidType = fread(fID, 1, '*uint8');
fseek(fID, 5, 'cof');
dimensionCount = fread(fID, 1, '*uint32');
for ii = 1 : dimensionCount
    dimension = fread(fID, 4, '*char');
    sBlockHeader.([dimension(1),'Start']) = fread(fID, 1, '*uint32');
    if ~strcmp(dimension(1),'X') && ~strcmp(dimension(1),'Y')
        fseek(fID, 12, 'cof');
    else
        sBlockHeader.([dimension(1),'Size']) = fread(fID, 1, '*uint32');
        fseek(fID, 8, 'cof');
    end
end
sBlockHeader.dataStartPos = currPos + 256 + mDataSize;
end

function fileHeader = readFILEHeader(fID)
fileHeader.major = fread(fID, 1, '*uint32');
fileHeader.minor = fread(fID, 1, '*uint32');
fseek(fID, 8, 'cof');
fileHeader.primaryFileGuid = fread(fID, 2, '*uint64');
fileHeader.fileGuid = fread(fID, 2, '*uint64');
fileHeader.filePart = fread(fID, 1, '*uint32');
fileHeader.dirPos = fread(fID, 1, '*uint64');
fileHeader.mDataPos = fread(fID, 1, '*uint64');
fseek(fID, 4, 'cof');
fileHeader.attDirPos  = fread(fID, 1, '*uint64');
end

function readAttach(fID)
dataSize = fread(fID, 1, '*uint32');
fseek(fID, 24, 'cof');
filePos = fread(fID, 1, '*uint64');
fseek(fID, 20, 'cof');
contentType = fread(fID, 8, '*char')';
disp(contentType)
name = fread(fID, 80, '*char')';
disp(name)
if strfind(contentType, 'JPG')
    fseek(fID, 112, 'cof');
    fout = fopen('thumbnail.jpg', 'wb');
    fwrite(fout, fread(fID, dataSize, '*uint8'), 'uint8');
    fclose(fout);
end
end

function pixType = getPixType(index)
switch index
    case 0
        pixType = 'Gray8';
    case 1
        pixType = 'Gray16';
    case 2
        pixType = 'Gray32Float';
    case 3
        pixType = 'Bgr24';
    case 4
        pixType = 'Bgr48';
    case 8
        pixType = 'Bgr96Float';
    case 9
        pixType = 'Bgra32';
    case 10
        pixType = 'Gray64ComplexFloat';
    case 11
        pixType = 'Bgr192ComplexFloat';
    case 12
        pixType = 'Gray32';
    case 13
        pixType = 'Gray64';
end
end

function compType = getCompType(index)
if index >= 1000
    compType = 'System-RAW';
elseif index >= 100 && index < 999
    compType = 'Camera-RAW';
else 
    switch index
        case 0
            compType = 'Uncompressed';
        case 1
            compType = 'JPEG';
        case 2
            compType = 'LZW';
        case 4
            compType = 'JPEG-XR';
    end
end
end

%% imshow 3D

function  imshow3D( Img, disprange, initS )
%IMSHOW3D displays 3D grayscale or RGB images in a slice by slice fashion
%with mouse-based slice browsing and window and level adjustment control,
%and auto slice browsing control.
%
% Usage:
% imshow3D ( Image )
% imshow3D ( Image , [] )
% imshow3D ( Image , [LOW HIGH] )
% imshow3D ( Image , [] , initsn )
%   
%    Image:      3D image MxNxKxC (K slices of MxN images) C is either 1
%                (for grayscale images) or 3 (for RGB images)  
%    [LOW HIGH]: display range that controls the display intensity range of
%                a grayscale image (default: the broadest available range)
%    initsn:     The slice number to be displayed initially (default:
%                mid-slice number) 
%
% Use the scroll bar or mouse scroll wheel to switch between slices. To
% adjust window and level values keep the mouse right button pressed, and
% drag the mouse up and down (for level adjustment) or right and left (for
% window adjustment). Window and level adjustment control works only for
% grayscale images.
% "Play" button displays all the slices as a sequence of frames. The time
% interval value can also be adjusted (default time interval is 100 ms).
% 
% "Auto W/L" button adjust the window and level automatically for grayscale
% images.
%
% While "Fine Tune" checkbox is checked the window/level adjustment gets 16
% times less sensitive to mouse movement, to make it easier to control
% display intensity rang.
%
% Note: The sensitivity of mouse-based window and level adjustment is set
% based on the user-defined display intensity range; the wider the range,
% the more sensitivity to mouse drag.
% 
% Note: IMSHOW3DFULL is a newer version of IMSHOW3D (also available on
% MathWorks) that displays 3D grayscale or RGB images from three
% perpendicular views (i.e., axial, sagittal, and coronal).
% 
%   Example
%   --------
%       % To display an image (MRI example)
%       load mri 
%       Image = squeeze(D); 
%       figure, 
%       imshow3D(Image) 
%
%       % To display the image, and adjust the display range
%       figure,
%       imshow3D(Image,[20 100]);
%
%       % To define the initial slice number
%       figure,
%       imshow3D(Image,[],5);
%
%   See also IMSHOW.
%
% - Maysam Shahedi (mshahedi@gmail.com)
% - Released: 1.0.0   Date: 2013/04/15
% - Revision: 1.1.0   Date: 2013/04/19
% - Revision: 1.5.0   Date: 2016/09/22
% - Revision: 1.6.0   Date: 2018/06/07
% - Revision: 1.6.1   Date: 2018/10/29
% 
sno = size(Img,3);  % number of slices
S = round(sno/2);
PlayFlag = false;   % Play flag, playing when it is 'True'
Tinterv = 100;
global InitialCoord;
MinV = 0;
MaxV = max(Img(:));
LevV = (double( MaxV) + double(MinV)) / 2;
Win = double(MaxV) - double(MinV);
WLAdjCoe = (Win + 1)/1024;
FineTuneC = [1 1/16];    % Regular/Fine-tune mode coefficients
if isa(Img,'uint8')
    MaxV = uint8(Inf);
    MinV = uint8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint16')
    MaxV = uint16(Inf);
    MinV = uint16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint32')
    MaxV = uint32(Inf);
    MinV = uint32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint64')
    MaxV = uint64(Inf);
    MinV = uint64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int8')
    MaxV = int8(Inf);
    MinV = int8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int16')
    MaxV = int16(Inf);
    MinV = int16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int32')
    MaxV = int32(Inf);
    MinV = int32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int64')
    MaxV = int64(Inf);
    MinV = int64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'logical')
    MaxV = 0;
    MinV = 1;
    LevV =0.5;
    Win = 1;
    WLAdjCoe = 0.1;
end    
SFntSz = 9;
txtFntSz = 10;
LVFntSz = 9;
WVFntSz = 9;
BtnSz = 10;
if (nargin < 3)
    S = round(sno/2);
else
    S = initS;
    if S > sno
        S = sno;
        warning('Initial slice number out of range');
    elseif S < 1
        S = 1;
        warning('Initial slice number out of range');
    end
end
if (nargin < 2)
    [Rmin Rmax] = WL2R(Win, LevV);
elseif numel(disprange) == 0
    [Rmin Rmax] = WL2R(Win, LevV);
else
    LevV = (double(disprange(2)) + double(disprange(1))) / 2;
    Win = double(disprange(2)) - double(disprange(1));
    WLAdjCoe = (Win + 1)/1024;
    [Rmin Rmax] = WL2R(Win, LevV);
end
clf
axes('position',[0,0.2,1,0.8]), imshow(squeeze(Img(:,:,S,:)), [Rmin Rmax])
FigPos = get(gcf,'Position');
S_Pos = [30 45 uint16(FigPos(3)-100)+1 20];
Stxt_Pos = [30 65 uint16(FigPos(3)-100)+1 15];
Wtxt_Pos = [20 18 60 20];
Wval_Pos = [75 20 50 20];
Ltxt_Pos = [130 18 45 20];
Lval_Pos = [170 20 50 20];
Btn_Pos = [240 20 70 20];
ChBx_Pos = [320 20 80 20];
Play_Pos = [uint16(FigPos(3)-100)+40 45 30 20];
Time_Pos = [uint16(FigPos(3)-100)+35 20 40 20];
Ttxt_Pos = [uint16(FigPos(3)-100)-50 18 90 20];
% W/L Button styles:
WL_BG = ones(Btn_Pos(4),Btn_Pos(3),3)*0.85;
WL_BG(1,:,:) = 1; WL_BG(:,1,:) = 1; WL_BG(:,end-1,:) = 0.6; WL_BG(:,end,:) = 0.4; WL_BG(end,:,:) = 0.4;
% Play Button styles:
Play_BG = ones(Play_Pos(4),Play_Pos(3),3)*0.85;
Play_BG(1,:,:) = 1; Play_BG(:,1,:) = 1; Play_BG(:,end-1,:) = 0.6; Play_BG(:,end,:) = 0.4; Play_BG(end,:,:) = 0.4;
Play_Symb = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1;...
             0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1;...
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1;...
             0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;...
             0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
Play_BG(floor((Play_Pos(4)-13)/2)+1:floor((Play_Pos(4)-13)/2)+13,floor(Play_Pos(3)/2)-7:floor(Play_Pos(3)/2)+6,:) = ...
    repmat(Play_Symb,1,1,3) .* Play_BG(floor((Play_Pos(4)-13)/2)+1:floor((Play_Pos(4)-13)/2)+13,floor(Play_Pos(3)/2)-7:floor(Play_Pos(3)/2)+6,:);
Pause_BG = ones(Play_Pos(4),Play_Pos(3),3)*0.85;
Pause_BG(1,:,:) = 1; Pause_BG(:,1,:) = 1; Pause_BG(:,end-1,:) = 0.6; Pause_BG(:,end,:) = 0.4; Pause_BG(end,:,:) = 0.4;
Pause_Symb = repmat([0, 0, 0, 1, 1, 1, 1, 0, 0, 0],13,1);
Pause_BG(floor((Play_Pos(4)-13)/2)+1:floor((Play_Pos(4)-13)/2)+13,floor(Play_Pos(3)/2)-5:floor(Play_Pos(3)/2)+4,:) = ...
    repmat(Pause_Symb,1,1,3) .* Pause_BG(floor((Play_Pos(4)-13)/2)+1:floor((Play_Pos(4)-13)/2)+13,floor(Play_Pos(3)/2)-5:floor(Play_Pos(3)/2)+4,:);
if sno > 1
    shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'FontSize', SFntSz);
    playhand = uicontrol('Style', 'pushbutton','Position', Play_Pos, 'Callback' , @Play);
    set(playhand, 'cdata', Play_BG)
    ttxthand = uicontrol('Style', 'text','Position', Ttxt_Pos,'String','Interval (ms): ',  'FontSize', txtFntSz);
    timehand = uicontrol('Style', 'edit','Position', Time_Pos,'String',sprintf('%d',Tinterv), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @TimeChanged);
else
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'FontSize', SFntSz);
end    
ltxthand = uicontrol('Style', 'text','Position', Ltxt_Pos,'String','Level: ',  'FontSize', txtFntSz);
wtxthand = uicontrol('Style', 'text','Position', Wtxt_Pos,'String','Window: ',  'FontSize', txtFntSz);
lvalhand = uicontrol('Style', 'edit','Position', Lval_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @WinLevChanged);
wvalhand = uicontrol('Style', 'edit','Position', Wval_Pos,'String',sprintf('%6.0f',Win), 'BackgroundColor', [1 1 1], 'FontSize', WVFntSz,'Callback', @WinLevChanged);
Btnhand = uicontrol('Style', 'pushbutton','Position', Btn_Pos,'String','Auto W/L', 'FontSize', BtnSz, 'Callback' , @AutoAdjust);
set(Btnhand, 'cdata', WL_BG)
ChBxhand = uicontrol('Style', 'checkbox','Position', ChBx_Pos,'String','Fine-tune', 'FontSize', txtFntSz);
set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
set (gcf, 'ButtonDownFcn', @mouseClick);
set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
set(gcf,'WindowButtonUpFcn', @mouseRelease)
set(gcf,'ResizeFcn', @figureResized)
% -=< Figure resize callback function >=-
    function figureResized(object, eventdata)
        FigPos = get(gcf,'Position');
        S_Pos = [30 45 uint16(FigPos(3)-100)+1 20];
        Stxt_Pos = [30 65 uint16(FigPos(3)-100)+1 15];
        Play_Pos = [uint16(FigPos(3)-100)+40 45 30 20];
        Time_Pos = [uint16(FigPos(3)-100)+35 20 40 20];
        Ttxt_Pos = [uint16(FigPos(3)-100)-50 18 90 20];
        if sno > 1
            set(shand,'Position', S_Pos);
            set(playhand, 'Position', Play_Pos)
            set(ttxthand, 'Position', Ttxt_Pos)
            set(timehand, 'Position', Time_Pos)
        end
        set(stxthand,'Position', Stxt_Pos);
        set(ltxthand,'Position', Ltxt_Pos);
        set(wtxthand,'Position', Wtxt_Pos);
        set(lvalhand,'Position', Lval_Pos);
        set(wvalhand,'Position', Wval_Pos);
        set(Btnhand,'Position', Btn_Pos);
        set(ChBxhand,'Position', ChBx_Pos);
    end
% -=< Slice slider callback function >=-
    function SliceSlider (hObj,event, Img)
        S = round(get(hObj,'Value'));
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
    end
% -=< Mouse scroll wheel callback function >=-
    function mouseScroll (object, eventdata)
        UPDN = eventdata.VerticalScrollCount;
        S = S - UPDN;
        if (S < 1)
            S = 1;
        elseif (S > sno)
            S = sno;
        end
        if sno > 1
            set(shand,'Value',S);
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
    end
% -=< Mouse button released callback function >=-
    function mouseRelease (object,eventdata)
        set(gcf, 'WindowButtonMotionFcn', '')
    end
% -=< Mouse click callback function >=-
    function mouseClick (object, eventdata)
        MouseStat = get(gcbf, 'SelectionType');
        if (MouseStat(1) == 'a')        %   RIGHT CLICK
            InitialCoord = get(0,'PointerLocation');
            set(gcf, 'WindowButtonMotionFcn', @WinLevAdj);
        end
    end
% -=< Window and level mouse adjustment >=-
    function WinLevAdj(varargin)
        PosDiff = get(0,'PointerLocation') - InitialCoord;
        Win = Win + PosDiff(1) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        LevV = LevV - PosDiff(2) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        if (Win < 1)
            Win = 1;
        end
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
        InitialCoord = get(0,'PointerLocation');
    end
% -=< Window and level text adjustment >=-
    function WinLevChanged(varargin)
        LevV = str2double(get(lvalhand, 'string'));
        Win = str2double(get(wvalhand, 'string'));
        if (Win < 1)
            Win = 1;
        end
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
    end
% -=< Window and level to range conversion >=-
    function [Rmn Rmx] = WL2R(W,L)
        Rmn = L - (W/2);
        Rmx = L + (W/2);
        if (Rmn >= Rmx)
            Rmx = Rmn + 1;
        end
    end
% -=< Window and level auto adjustment callback function >=-
    function AutoAdjust(object,eventdata)
        Win = double(max(Img(:))-min(Img(:)));
        Win (Win < 1) = 1;
        LevV = double(min(Img(:)) + (Win/2));
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
    end
% -=< Play button callback function >=-
    function Play (hObj,event)
        PlayFlag = ~PlayFlag;
        if PlayFlag
            set(playhand, 'cdata', Pause_BG)
        else
            set(playhand, 'cdata', Play_BG)
        end            
        while PlayFlag
            S = S + 1;
            if (S > sno)
                S = 1;
            end
            set(shand,'Value',S);
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
            set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
            pause(Tinterv/1000)
        end
    end
% -=< Time interval adjustment callback function>=-
    function TimeChanged(varargin)
        Tinterv = str2double(get(timehand, 'string'));
    end
    
end
% -=< Maysam Shahedi (mshahedi@gmail.com), October 29, 2018>=-
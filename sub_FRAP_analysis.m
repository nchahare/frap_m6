%% Created on
%   09-Aug-2019 10:28:51
%% @author: nchahare
clear all; clc; close all

%% <<<<<<<<<<<<<<<<<<< description
% this is for analysing the frap data. To browse the script, use "Go To"

% load('bigdata.mat')

[file,path] = uigetfile('*.mat');
load(fullfile(path,file));

for imagenumber = 1:length(RAWDATA)

%% load data in the basic parameters

imgs = RAWDATA(imagenumber).Images;
% mask = double(RAWDATA(imagenumber).Mask);

medianI = double(extractfield(RAWDATA(imagenumber).FrapIntensity,'median'));
time = extractfield(RAWDATA(imagenumber).FrapIntensity,'time');

Rscale = sqrt(2)*RAWDATA(imagenumber).xscale;


%% find the frap image

frapped_pt = 21;
fraploc = RAWDATA(imagenumber).FrapCenter;
R = RAWDATA(imagenumber).FrapRadius;

justbeforefrapimg = mean(imgs(:,:,frapped_pt-3:frapped_pt-1),3);
frapedimg = double(imgs(:,:,frapped_pt));


%% normalization of the time response

FI = double(extractfield(RAWDATA(imagenumber).FrapIntensity,'IntDen'));
CI = double(extractfield(RAWDATA(imagenumber).CellIntensity,'IntDen'));
time  = extractfield(RAWDATA(imagenumber).CellIntensity,'time');

F.rawprebleach = median(FI(1:frapped_pt-1));

% figure; plot(time,FI)

F.raw = doublenormalize(FI,CI,frapped_pt);

% figure; plot(time,F.raw,'b')

F.prebleach = median(F.raw(1:frapped_pt-1));

F.raw = F.raw(frapped_pt:end);
time  = time(frapped_pt:end);

%% filtering the response of the time response

% figure; plot(time,F.raw,'b')

% with movmean
framelen = 4;
F.MovingMean = movmean(F.raw,framelen);

% with Satvitzky Golay filter
order = 1; framelen = 5;
F.SatvitzkyGolay = sgolayfilt(F.raw,order,framelen);
hold on
plot(time,F.MovingMean,'g')

hold on
plot(time,F.SatvitzkyGolay,'r')


%% Quality check for the data


quality = 'good';
if bleachingdepth(FI,frapped_pt)<0.6
    quality = 'bad ';
end
if mobilefraction(F.raw) < 0
    quality = 'bad ';
end
if mobilefraction(F.MovingMean) < 0
    quality = 'bad ';
end
if mobilefraction(F.SatvitzkyGolay) < 0
    quality = 'bad ';
end

%% saving data into a structure 


ResutingData(imagenumber).file              = RAWDATA(imagenumber).File;
ResutingData(imagenumber).FrapRadius        = RAWDATA(imagenumber).FrapRadius;
ResutingData(imagenumber).mobileFractionSG  = abs(mobilefraction(F.SatvitzkyGolay));
ResutingData(imagenumber).mobileFractionMA  = abs(mobilefraction(F.MovingMean));
ResutingData(imagenumber).mobileFractionRaw = abs(mobilefraction(F.raw));
ResutingData(imagenumber).BD                = bleachingdepth(FI,frapped_pt);
ResutingData(imagenumber).Quality           = quality;



% figure
% subplot(1,2,1)
% plot(FI)
% subplot(1,2,2)
% plot(time,F.raw,'b')
% hold on
% plot(time,F.MovingMean,'g')
% hold on
% plot(time,F.SatvitzkyGolay,'r')
% title([num2str(bleachingdepth(FI,frapped_pt)) '---'   num2str(mobilefraction(F.SatvitzkyGolay))])

end

% save file as user input name
prompt = {'Save data file as:'};
dlgtitle = 'Input';
dims = [1 35];
filename = [char(inputdlg(prompt,dlgtitle,dims)) '.mat'];
save(filename,'RAWDATA')

%% --------------------------FUNCTIONS-------------------------------------
% normalization of the time response
function ynormal = doublenormalize(ydata,y0data, fraptime)

ypre  = median(ydata(1:fraptime-1));
y0pre = median(y0data(1:fraptime-1));

ynormal = ydata./ypre;
ynormal = ynormal.*(y0pre./y0data);

end

%calculate the mobile fraction
function mf = mobilefraction(ydata)

mf = (median(ydata(end-20:end)) - ydata(1))/(1 - ydata(1));

end

%calculate the bleaching depth
function bd = bleachingdepth(ydata,f)

bd = (ydata(f-1)- ydata(f))/ydata(f-1);

end
% Test out the reading and processesing of a Dicom File

filename = '6x_10x10_100SSD_Dose_PLAN.dcm';

info = dicominfo(filename);
I = dicomread(filename);

rows = double(info.Rows);
cols = double(info.Columns);
deps = double(info.NumberOfFrames);

% Rows Appear to be coronal
% Columns Appear to be Sagittal

x = info.PixelSpacing(1)*(((1:rows) - rows/2) - 0.5) ;
y = info.PixelSpacing(2)*(((1:cols) - cols/2) - 0.5) ;
z = info.PixelSpacing(2)*(((1:deps) - deps/2) - 0.5) ;

J = reshape(I(:,:,1,:),[rows cols deps]);

DOSE = double(J)*info.DoseGridScaling;

% I want a plane at: x = 0, y = 0 and z = 0
xloc = 0;
yloc = 40;
zloc = 0;

figure(1)


subplot(1,3,1);
DOSE2D = getPlaneAt(x,y,z,DOSE,xloc,'x');
axis('square')
imagesc(y,z,DOSE2D)

subplot(1,3,2);
DOSE2D = getPlaneAt(x,y,z,DOSE,yloc,'y');
imagesc(x,z,DOSE2D)

subplot(1,3,3);
DOSE2D = getPlaneAt(x,y,z,DOSE,zloc,'z');
imagesc(z,y,DOSE2D)

%% Test out the reading and processing of Acsii files from OmniPro

% Read one file from OmniPro and Create Structure
filename = 'P06_Open_OPD_2.ASC';
omniproStruct = omniproAccessTOmat(filename);

% Read another file and add to previous structure
filename = 'P06_Open_OPP.ASC';
omniproStruct = omniproAccessTOmat(filename,omniproStruct);

figure(2)
hold on;
% Get OPD (Open field PDD) for 15x15 field (150 mm by 150 mm)
[ x, y, z, d ] = getOmniproAccessData(omniproStruct,'OPD', [150 150]);
plot(z,d)
% Get OPP (Open field profile) crossline (X) profile for 15x15 field (150
% mm by 150 mm) at a depth of 5 cm (50 mm)
[ x, y, z, d ] = getOmniproAccessData(omniproStruct,'OPP', [150 150],50,'X');
plot(x,d)


% Read Jeni's file from OmniPro and Create Structure
filename = 'P06OPN_WISC.ASC';
omniproStruct2 = omniproAccessTOmat(filename);

% Get OPD (Open field PDD) for 10x10 field (100 mm by 100 mm)
[ x, y, z, d ] = getOmniproAccessData(omniproStruct2,'OPD', [100 100]);
plot(z,d,'r')
hold off;

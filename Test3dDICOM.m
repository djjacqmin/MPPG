filename = 'SampleDose.dcm';

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
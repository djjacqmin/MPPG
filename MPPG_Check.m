function MPPG_Check
%MPPG_Check.m
%This script compares measured and calculated MPPG data for all tests and energies.
clear all;
close all;

meas_f_path = 'E:\\Google Drive\\Medical_Physics\\MPPG\\MPPG_BasicPhotonFields-2014-01-17\\';
meas_f_file = '20131223_MPPG_Meas_Data.xlsx';
meas_f_full = [meas_f_path meas_f_file];
[~,~,measTable] = xlsread(meas_f_full, 5);

mD = measTable(2:end,:); %remove header row
[numTests numCols] = size(mD);

plotOn = 1;

%Which columns hold what data in table, this makes it a little easier to change which columns are which
mPathID = 10;
mFnameID = 11;
mPageID = 4;
mDepRID = 5;
mIndepRID = 6;
testID = 1;
enID = 2;
subTestID = 3;
cPathID = 8;
cFnameID = 9;
cOrigXID = 12;
cOrigYID = 13;
cOrigZID = 14;
cOrigResID = 15;

%Loop through all verification tests
for v = 1:numTests
    if strcmp('PDD',mD(v,3))
        %Run PDD verification
        measPDD = GetMeasPDD(mD{v,mPathID}, mD{v,mFnameID}, mD{v,mPageID},...
            mD{v,mDepRID}, mD{v,mIndepRID}, plotOn, ...
            mD{v,testID}, mD{v,enID}, mD{v,subTestID});
        
        calcPDD = GetCalcPDD(mD{v,cPathID}, mD{v,cFnameID}, mD{v,cOrigXID},...
            mD{v,cOrigYID}, mD{v,cOrigZID}, mD{v,cOrigResID}, plotOn, ...
            mD{v,testID}, mD{v,enID}, mD{v,subTestID});
        %[regMeas regCalc] = RegisterPDDs(measPDD, calcPDD);
        %vOut = VerifyPDD(regMeas, regCalc);
        %WriteSummary(vOut);
    else
        %Run cross beam verification
        
    end
        
end

function measPDD = GetMeasPDD(path, fname, page, depR, indepR, plotOn, test, en, sub)
    dep = xlsread([path fname], page, depR); %dependent var
    indep = xlsread([path fname], page, indepR);
    measPDD = [indep dep];
    
    if plotOn
        figure;
        plot(indep,dep);
        ttl = ['Measured ' test ' ' num2str(en) 'MV ' sub];
        title(ttl);
    end
end

function calcPDD = GetCalcPDD(path, fname, origX, origY, origZ, res, plotOn, test, en, sub)
    
    info = dicominfo([path fname]);    
    Y = dicomread(info);
    X = squeeze(Y(:,:,1,:));
    
    %Get size of array
    [numPixX numPixY numPixZ] = size(X);

    %Get interogation point in image frame
    %The image frame has 0,0,0 at most right, post, sup
    pixX = round((-origX)/res);
    pixSurY = 1; %water surface in image frame
    pixZ = numPixZ - round(-origZ/res);    
    OnAxPdd = squeeze(squeeze(X(pixSurY:end,pixX,pixZ)));
    OnAxPdd = double(OnAxPdd)/max(double(OnAxPdd));
    PddIndep = linspace(0,length(OnAxPdd)*res,length(OnAxPdd));
    %interpolate
    interpFac = 10;
    PddIndepInt = linspace(0,length(OnAxPdd)*res,length(OnAxPdd)*interpFac);
    OnAxPddInt = spline(PddIndep, OnAxPdd, PddIndepInt);

    if plotOn
        figure;
        plot(PddIndep,OnAxPdd); hold all;
        plot(PddIndepInt,OnAxPddInt);
        ttl = ['Calculated ' test ' ' num2str(en) 'MV ' sub];
        title(ttl);        
    end
    
    calcPDD = 0;
end

function [regMeas regCalc] = RegisterPDDs(measPDD, calcPDD)
end

end
%Modification notes (CYF)
%111414 limited T1 and M0 conversion to standard DICOMs to 1st image
%022515 output GMprob, meanEPI and T1 in CBF_folder for coregistration with
%lesion masks
%082615 added dicom header info to settings so they can be passed to markdown file
function [gGMCBF,SNRout,outliers]=asl_process100317(aslfolder,t1folder,M0folder,Lesionmask,age,fieldstrength,subtrmethod,outliermethod,kernel,fcnl,aslprotocol,w,tau,acq,NLM,PVcmethod)
warning off all;close all;
global defaults;
spm_defaults;
global settings;

setenv( 'FSLDIR', '/home/tbp688/fmri/fsl');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;

sprintf('Processing folder: %s\n',aslfolder)
    cd(aslfolder)
    settings.CBF_folder=fullfile(aslfolder,'CBFmaps');
    if ~exist('settings.CBF_folder','dir')
        status=mkdir(settings.CBF_folder);
    end
    
if isempty(t1folder)
    PVcmethod='none';
end

settings.field=fieldstrength;%input('Enter field strength: 1.5, 3 or 7\n');
switch settings.field
    case 1.5, settings.T1b=1300;
    case 3, settings.T1b=1664;
    case 7, settings.T1b=2200;
end
settings.regress=1;
settings.subt=subtrmethod;%input('Select subtraction method: 1=pairwise, 2=surround\n');
%settings.lambda=0.9;%g/ml
if kernel~=0
    settings.smooth=1;settings.kernel=kernel;
else
    settings.smooth=0;
end

settings.FSL=0;

settings.outliermethod=outliermethod;

settings.reslcflags=struct('interp',1,'mask',1,'mean',1,'which',2,'wrap',[0 0 0]');
%settings.fcnl=fcnl;%input('Output difference timeseries for ASL fMRI? 0=No, 1=yes\n');
%settings.removeoutlier=1;
if (strcmp(aslprotocol,'PASL'))
    settings.asl=1;settings.eff=0.98;else
    settings.asl=2;settings.eff=0.85;end
settings.tau=tau;settings.w=w;

settings.WIP=0;settings.NLM=NLM;
settings.acq=acq;
settings.PVcmethod=PVcmethod;
settings.is3D=0;
%extract info from dicoms, convert them to NIFTI for SPM scripts
try
    files=spm_select('FPList',aslfolder,'.*');
    
    doscmd=sprintf('gdcm gdcminfo %s',files(1,:));[status,cmdout]=dos(doscmd);
    %Convert enhanced dicom to standard dicom
    if strfind(cmdout,'Enhanced')
        display('Enhanced DICOM detected. Converting to standard DICOMS...\n');
        fullfolder=fullfile(aslfolder,'Standard');
        if exist(fullfolder,'dir')~=7
            mkdir(fullfolder);
        end
        
        doscmd=sprintf('emf2sf --out-dir %s %s',fullfolder,files(1,:));
        dos(doscmd);
        aslfolder=fullfolder;
    end
    
    %Get header info & convert to NIFTI
    cd(aslfolder);disp('In ASL folder...\n');
    files=spm_select('FPList',pwd,'.*');
    if size(files,1)>20
        settings.outlier=1;
    else settings.outlier=0;end
    display('Reading ASL headers...\n');
    hdr=spm_dicom_headers(files);
    display('Done!\n');
    settings.serNumber=hdr{1,1}.SeriesNumber;
    settings.ser = regexprep(hdr{1,1}.SeriesDescription,'\s','_');
    settings.Voxels=[hdr{1,1}.PixelSpacing;hdr{1,1}.SliceThickness];
    settings.TR=hdr{1,1}.RepetitionTime;
    settings.TE=hdr{1,1}.EchoTime;
    settings.ScanDate=datestr(hdr{1,1}.StudyDate,'mm/dd/yyyy');
    if isfield(hdr{1,1},'PatientsName')
        settings.PatientID=deblank(hdr{1,1}.PatientsName);
    else
        settings.PatientID=deblank(hdr{1,1}.PatientName);
    end
    if isfield(hdr{1,1},'PatientsAge')
        settings.Age=hdr{1,1}.PatientsAge;
    elseif isfield(hdr{1,1},'PatientAge')
        settings.Age=hdr{1,1}.PatientAge
    else
        settings.Age=age;
    end
    if isfield(hdr{1,1},'PatientsSex')
        settings.gender=hdr{1,1}.PatientsSex;
    else
        settings.gender=hdr{1,1}.PatientSex;
    end
    if strcmp(settings.PatientID,'')
        settings.PatientID='ID';
    end
    settings.visit=deblank(hdr{1,1}.PatientID);
    settings.Site=hdr{1,1}.InstitutionName;
    spm_dicom_convertP50(hdr,'all','flat','nii');
    
    
    if strfind(settings.PatientID,'BU')%change to Acquisition Site?
        settings.w=2000;
    end
    
    
    if strcmp(hdr{1,1}.MRAcquisitionType,'3D')
        settings.acq=0;settings.regress=0;%settings.NLM=0;
        settings.is3D=1;
        %settings.eff=settings.eff*0.75;%BS adjustment based on Garcia 2005 paper
        if (~isempty(strfind(hdr{1,1}.Manufacturer,'SIEMENS'))&&~isempty(strfind(hdr{1,1}.SequenceName,'tg818')))
            settings.WIP=1;
        end
    end
    
    %CYF added 3/31/14-----------------------------------
    %extract motion parameters from Siemens moco series image comment
    if ~isempty(strfind(settings.ser,'MoCo'))
        tempmat=zeros(size(files,1)-1,6);mocomat=tempmat;
        for b=2:size(files,1)
            info=dicominfo(files(b,:));
            if strfind(info.ImageComments,'Motion')
                if isempty(strfind(info.ImageComments,'Pro'))
                    temp=info.ImageComments(9:end);
                else
                    temp=info.ImageComments(9:(strfind(info.ImageComments,'Pro')-1));
                end
                tempmat(b,:)=str2num(temp);end
        end
        niftis=spm_select('FPList',aslfolder,'^f.*');
        [pth,nm]=fileparts(niftis(1,:));
        fname=fullfile(pth,['rp_',nm,'.txt']);
        mocomat=[tempmat(:,2),tempmat(:,1),-tempmat(:,3),-tempmat(:,5)/(180/pi),tempmat(:,4)/(180/pi),tempmat(:,6)/(180/pi)];
        save(fname,'mocomat','-ascii');
    end
    
    
    %convert T1s
    if ~isempty(t1folder)
        cd(t1folder);disp('In T1 folder...\n')
        T1files=spm_select('FPList',pwd,'.*');
        [pth,nm,ext]=fileparts(T1files(1,:));                                  
        if strcmp(ext,'.nii')~=1 
        doscmd=sprintf('gdcm gdcminfo %s',deblank(T1files(1,:)));[status,cmdout]=dos(doscmd);
        %Convert enhanced dicom to standard dicom
        if strfind(cmdout,'Enhanced')
            disp('Enhanced DICOM detected. Converting to standard DICOMS...\n');
            fullfname=fullfile(t1folder,'Standard');
            if exist(fullfname,'dir')~=7
                mkdir(fullfname);
            end
            doscmd=sprintf('emf2sf --out-dir %s %s',fullfname,deblank(T1files(1,:)));
            dos(doscmd);
            t1folder=fullfname;
            cd(t1folder)
            T1files=spm_select('FPList',pwd,'.*');
        end
        hdr=spm_dicom_headers(T1files);
        spm_dicom_convertP50(hdr,'all','flat','nii');
        end
    end
    
    if ~isempty(M0folder)
        cd(M0folder); disp('In M0 folder...\n')
        M0files=spm_select('FPList',pwd,'.*');
        doscmd=sprintf('gdcm gdcminfo %s',deblank(M0files(1,:)));[status,cmdout]=dos(doscmd);
        %Convert enhanced dicom to standard dicom
        if strfind(cmdout,'Enhanced')
            display('Enhanced DICOM detected. Converting to standard DICOMS...\n');
            fullfname=fullfile(M0folder,'Standard');
            if exist(fullfname,'dir')~=7
                mkdir(fullfname);
            end
            doscmd=sprintf('emf2sf --out-dir %s %s',fullfname,deblank(M0files(1,:)));
            dos(doscmd);
            M0folder=fullfname;
            cd(M0folder)
            M0files=spm_select('FPList',pwd,'.*');
        end
        hdr=spm_dicom_headers(M0files);
        spm_dicom_convertP50(hdr,'all','flat','nii');
        
    end
    
    spm_figure('CreateWin','Interactive');
    F=spm_figure('FindWin','Interactive');set(F,'Position',[50 50 500 500]);
    spm_figure('CreateWin','Graphics');
    settings.ext='nii';
    
    
    settings.pdffolder=fullfile(settings.CBF_folder,'PDFs');
    if ~exist(settings.pdffolder,'dir')
        status=mkdir(settings.pdffolder);
    end
    
    settings.intfolder=fullfile(settings.CBF_folder,'Intermediates');
    if ~exist(settings.intfolder,'dir')
        status=mkdir(settings.intfolder);
    end
    
    [gGMCBF,SNRout,outliers,f,Slices,isquant]=asl_pipeline100317(aslfolder,t1folder,M0folder,Lesionmask,fcnl,settings,defaults);
    %     if settings.FSL==0
    %         [gGMCBF]=asl_quantification092217(Perf,M0,mask,V_CBF,settings,defaults);
    %     else
    %         %TO DO: input FSL commands
    %         %TO DO: read in FSL CBF results
    %     end
    
    %Kate: Please update location of this file.
    %fID=fopen('E:\NUdata\Swinburne/summary.txt','a+')
    fID=fopen(fullfile(aslfolder,'params.ini'),'w')
    fprintf(fID,'subtraction method=%s\n', settings.subt);
    if settings.smooth==1
        fprintf(fID,'smoothing kernel=%d\n', settings.kernel);
    else
        fprintf(fID,'smoothing kernel=0\n');
    end
    fprintf(fID,'rawPerfSNR=%4.2f\n', SNRout.rawPerfSNR);
    fprintf(fID,'temporalPerfSNR=%4.2f\n', SNRout.temporalPerfSNR);
    fprintf(fID,'rawCtrlSNR=%4.2f\n', SNRout.rawCtrlSNR);%added 070215
    fprintf(fID,'mCtrlSignal=%4.2f\n', SNRout.mCtrlSignal);%added 070215
    if length(gGMCBF)>1
        fprintf(fID,'GlobalGMCBF=%4.2f\n', gGMCBF(2));
    else
        fprintf(fID,'GlobalGMCBF=%4.2f\n', gGMCBF(1));
    end
    fprintf(fID,'deleteN=%d\n', outliers.dInd);
    fprintf(fID,'deletedPairs=%s\n', mat2str(outliers.DeletedPairs));
    fclose(fID)
    
    %% create Markdown file
    
    fileID=fopen(fullfile(aslfolder,['Report1.md']),'w');
    fprintf(fileID,'Report generated by ASLpipeline version %s in Matlab %s on %s\n',func2str(f),version, datestr(now));
    fprintf(fileID,'Please send questions & feedback to <a href=mailto:yfchen@northwestern.edu>Jennie Chen</a>.\n');
    
    fprintf(fileID,'\n111414 changed code to accomodate no T1 such that will run to completion.\n');
    fprintf(fileID,'070915 added option for regression of motion, 1st derivative of motion, gsignal & CSF.\n');
    fprintf(fileID,'082415 define mask using mean of ASL ctrls instead of mean, also implemented markdown for pdf generation.\n');
    fprintf(fileID,'031416 added functionality for Siemens WIP 818F BS-ASL. Modified qCBF code to remove NaN box around volume. Output indices of rejected volumes in param.ini.\n');
    fprintf(fileID,'100317 Rewrote pipeline to include denoising, PVc and normalization options.\n');
    fprintf(fileID,'121117 added option for LPCA DWIdenoising (Pierrick Coupe & Jose Manjon) of raw image time series.\n');
    
    fprintf(fileID,'\n##References\n');
    fprintf(fileID,'1. Alsop D., et al, MRM 2015;73:102-116\n');
    fprintf(fileID,'2. Tan H., et al, JMRI 2009;29:1134-1139\n ');
    fprintf(fileID,'3. Maumet C. et al., In: Multimodal Brain Image Analysis. MBIA 2012. Lecture Notes in Computer Science, vol 7509, pp215-224.\n ');
    fprintf(fileID,'4. Du AT., et al, Neurology 2006;67:1215-1220\n ');
    fprintf(fileID,'5. Asllani I., et al, MRM 2008;60:1362-1371\n ');
    fprintf(fileID,'6. Manjon J., et al, PLoS ONE 8(9): e73021\n ');
    fprintf(fileID,'7. Wang Z., et al, MRI 2012;30:1409-1415\n ');
    
    fprintf(fileID,'\n## ASLQA summary:\n');
    fprintf(fileID,'Name | Value\n');
    fprintf(fileID,':---------- | ----------:\n');
    if length(gGMCBF)>1
        fprintf(fileID,'Global GM CBF | %.2f ml/100g/min\n',gGMCBF(2));
    else
        fprintf(fileID,'Global GM CBF | %.2f ml/100g/min\n',gGMCBF(1));
    end
    fprintf(fileID,'rawPerfSNR | %.2f\n',SNRout.rawPerfSNR);
    fprintf(fileID,'tempPerfSNR | %.2f\n',SNRout.temporalPerfSNR);
    fprintf(fileID,'rawCtrlSNR | %.2f\n',SNRout.rawCtrlSNR);
    fprintf(fileID,'Ndeleted/Total | %d/%d\n',outliers.dInd,outliers.dInd+outliers.FinalPairNo);
    fclose(fileID);doscmd=sprintf('md2pdf Report1.md %s/Report1.pdf',settings.pdffolder);dos(doscmd)
    
    
    fileID=fopen(fullfile(aslfolder,['Report2.md']),'w');
    fprintf(fileID,'\n## ASLQA input parameters:\n');
    fprintf(fileID,'Name | Value\n');
    fprintf(fileID,':---------- | ----------:\n');
    fprintf(fileID,'FieldStrength | %d\n',settings.field);
    fprintf(fileID,'T1b | %d ms\n',settings.T1b);
    fprintf(fileID,'subtraction method | %s\n',settings.subt);
    fprintf(fileID,'Presmoothing | %d\n',settings.smooth);
    if settings.smooth~=0
        fprintf(fileID,'Presmoothing kernel | %d mm\n',settings.kernel);end
    fprintf(fileID,'fcnl | %d\n',fcnl);
    fprintf(fileID,'Outlier Removal | %s\n',settings.outliermethod);
    fprintf(fileID,'asl | %d\n',settings.asl);
    fprintf(fileID,'tau | %d\n',settings.tau);
    fprintf(fileID,'w | %d\n',settings.w);
    fprintf(fileID,'eff | %.2f\n',settings.eff);
    fprintf(fileID,'acq | %d\n',settings.acq);
    fprintf(fileID,'WIP | %d\n',settings.WIP);
    fprintf(fileID,'Regress motion&gsig | %d\n',settings.regress);
    if settings.NLM==1
        fprintf(fileID,'denoise | on\n');
    else
        fprintf(fileID,'denoise | off\n');
    end
        
    fclose(fileID);doscmd=sprintf('md2pdf Report2.md %s/Report2.pdf',settings.pdffolder);dos(doscmd)
    
    
    fileID=fopen(fullfile(aslfolder,['Report4.md']),'w');
    fprintf(fileID,'\n## Acquisition Parameters:\n');
    fprintf(fileID,'Name | Value\n');
    fprintf(fileID,':---------- | ----------:\n');
    fprintf(fileID,'SubjectID | %s\n',settings.PatientID);
    fprintf(fileID,'SubjectAge | %s\n',settings.Age);
    fprintf(fileID,'SubjectSex | %s\n',settings.gender);
    fprintf(fileID,'Visit | %s\n',settings.visit);
    fprintf(fileID,'Series | %03d_%s\n',settings.serNumber,settings.ser);
    fprintf(fileID,'SiteID | %s\n',settings.Site);
    fprintf(fileID,'Scan Date | %s\n',settings.ScanDate);
    fprintf(fileID,'Voxel Size | %.2f %.2f %.2f\n',settings.Voxels');
    fprintf(fileID,'Slices | %d\n',Slices);
    fprintf(fileID,'TR | %d\n',settings.TR);
    fprintf(fileID,'TE | %d\n',settings.TE);
    fclose(fileID);doscmd=sprintf('md2pdf Report4.md %s/Report4.pdf',settings.pdffolder);dos(doscmd)
    
    fileID=fopen(fullfile(aslfolder,['Report5.md']),'w');
    fprintf(fileID,'![SNR](SNR.png)\n');
    fprintf(fileID,'\n## SNR summary:\n');
    fprintf(fileID,'Name | Value\n');
    fprintf(fileID,':---------- | ----------:\n');
    fprintf(fileID,'mPerfSignal | %0.2g\n',SNRout.mPerfSignal);
    fprintf(fileID,'rawPerfNoise | %0.2g\n',SNRout.rawPerfNoise);
    fprintf(fileID,'uncorr SNR | %0.2g\n',SNRout.uncorrSNR);
    fprintf(fileID,'uncorr tSNR | %0.2g\n',SNRout.uncorrtSNR);
    fprintf(fileID,'rawPerfSNR | %0.2g\n',SNRout.rawPerfSNR);
    fprintf(fileID,'temporalPerfNoise | %0.2g\n',SNRout.temporalPerfNoise);
    fprintf(fileID,'temporalPerfSNR | %0.2g\n',SNRout.temporalPerfSNR);
    fprintf(fileID,'mCtrlSignal | %0.2g\n',SNRout.mCtrlSignal);
    fprintf(fileID,'rawCtrlNoise | %0.2g\n',SNRout.rawCtrlNoise);
    fprintf(fileID,'rawCtrlSNR | %0.2g\n',SNRout.rawCtrlSNR);
    fclose(fileID);doscmd=sprintf('md2pdf Report5.md %s/Report5.pdf',settings.pdffolder);dos(doscmd)
    
    reportf=fullfile(settings.CBF_folder,sprintf('%s_%04d_%s_ASLQA.pdf',settings.visit,settings.serNumber,settings.ser));
    
    doscmd=sprintf('pdfunite %s/Report1.pdf %s/Report2.pdf %s/Report4.pdf %s/Realign.pdf',settings.pdffolder,settings.pdffolder,settings.pdffolder,settings.pdffolder);
    
    if settings.outlier==1 && strcmp(settings.outliermethod,'Zthreshold')
        doscmd=sprintf('%s %s/Report3.pdf',doscmd,settings.pdffolder);
    end
    
    doscmd=sprintf('%s %s/Report5.pdf',doscmd,settings.pdffolder);
    
    if isquant==1
        doscmd=sprintf('%s %s/globalCBF.pdf',doscmd,settings.pdffolder);
        
        if strcmp(settings.PVcmethod,'None')~=1
            doscmd=sprintf('%s %s/globalCBFPVc.pdf',doscmd,settings.pdffolder);
        end
    else
        doscmd=sprintf('%s %s/meanPerf.pdf',doscmd,settings.pdffolder);
    end
    
    doscmd=sprintf('%s %s',doscmd,reportf);
    dos(doscmd)
    
    zip(fullfile(settings.CBF_folder,'intermediates.zip'),settings.intfolder);
    rmdir(settings.intfolder,'s')
    rmdir(settings.pdffolder,'s')
    
catch err
    disp(err.message);
    disp(err.identifier);
    for k=1:length(err.stack)
        fprintf('In %s at %d\n',err.stack(k).file,err.stack(k).line);
    end
    %exit;
end
return;




if (strcmp(aslprotocol,'PASL'))
    settings.asl=1;settings.eff=0.98;else
    settings.asl=2;settings.eff=0.85;end
settings.tau=tau;settings.w=w;

settings.WIP=0;settings.NLM=NLM;
settings.acq=acq;
settings.PVcmethod=PVcmethod;
settings.is3D=0;
%extract info from dicoms, convert them to NIFTI for SPM scripts
try
    files=spm_select('FPList',aslfolder,'.*');
    
    doscmd=sprintf('gdcm gdcminfo %s',files(1,:));[status,cmdout]=dos(doscmd);
    %Convert enhanced dicom to standard dicom
    if strfind(cmdout,'Enhanced')
        display('Enhanced DICOM detected. Converting to standard DICOMS...\n');
        fullfolder=fullfile(aslfolder,'Standard');
        if exist(fullfolder,'dir')~=7
            mkdir(fullfolder);
        end
        
        doscmd=sprintf('emf2sf --out-dir %s %s',fullfolder,files(1,:));
        dos(doscmd);
        aslfolder=fullfolder;
    end
    
    %Get header info & convert to NIFTI
    cd(aslfolder);disp('In ASL folder...\n');
    files=spm_select('FPList',pwd,'.*');
    if size(files,1)>20
        settings.outlier=1;
    else settings.outlier=0;end
    display('Reading ASL headers...\n');
    hdr=spm_dicom_headers(files);
    display('Done!\n');
    settings.serNumber=hdr{1,1}.SeriesNumber;
    settings.ser = regexprep(hdr{1,1}.SeriesDescription,'\s','_');
    settings.Voxels=[hdr{1,1}.PixelSpacing;hdr{1,1}.SliceThickness];
    settings.TR=hdr{1,1}.RepetitionTime;
    settings.TE=hdr{1,1}.EchoTime;
    settings.ScanDate=datestr(hdr{1,1}.StudyDate,'mm/dd/yyyy');
    if isfield(hdr{1,1},'PatientsName')
        settings.PatientID=deblank(hdr{1,1}.PatientsName);
    else
        settings.PatientID=deblank(hdr{1,1}.PatientName);
    end
    if isfield(hdr{1,1},'PatientsAge')
        settings.Age=hdr{1,1}.PatientsAge;
    elseif isfield(hdr{1,1},'PatientAge')
        settings.Age=hdr{1,1}.PatientAge
    else
        settings.Age=age;
    end
    if isfield(hdr{1,1},'PatientsSex')
        settings.gender=hdr{1,1}.PatientsSex;
    else
        settings.gender=hdr{1,1}.PatientSex;
    end
    if strcmp(settings.PatientID,'')
        settings.PatientID='ID';
    end
    settings.visit=deblank(hdr{1,1}.PatientID);
    settings.Site=hdr{1,1}.InstitutionName;
    spm_dicom_convertP50(hdr,'all','flat','nii');
    
    
    if strfind(settings.PatientID,'BU')%change to Acquisition Site?
        settings.w=2000;
    end
    
    
    if strcmp(hdr{1,1}.MRAcquisitionType,'3D')
        settings.acq=0;settings.regress=0;settings.NLM=0;settings.is3D=1;
        settings.eff=settings.eff*0.75;%BS adjustment based on Garcia 2005 paper
        if (~isempty(strfind(hdr{1,1}.Manufacturer,'SIEMENS'))&&~isempty(strfind(hdr{1,1}.SequenceName,'tg818')))
            settings.WIP=1;
        end
    end
    
    %CYF added 3/31/14-----------------------------------
    %extract motion parameters from Siemens moco series image comment
    if ~isempty(strfind(settings.ser,'MoCo'))
        tempmat=zeros(size(files,1)-1,6);mocomat=tempmat;
        for b=2:size(files,1)
            info=dicominfo(files(b,:));
            if strfind(info.ImageComments,'Motion')
                if isempty(strfind(info.ImageComments,'Pro'))
                    temp=info.ImageComments(9:end);
                else
                    temp=info.ImageComments(9:(strfind(info.ImageComments,'Pro')-1));
                end
                tempmat(b,:)=str2num(temp);end
        end
        niftis=spm_select('FPList',aslfolder,'^f.*');
        [pth,nm]=fileparts(niftis(1,:));
        fname=fullfile(pth,['rp_',nm,'.txt']);
        mocomat=[tempmat(:,2),tempmat(:,1),-tempmat(:,3),-tempmat(:,5)/(180/pi),tempmat(:,4)/(180/pi),tempmat(:,6)/(180/pi)];
        save(fname,'mocomat','-ascii');
    end
    
    
    %convert T1s
    if ~isempty(t1folder)
        cd(t1folder);disp('In T1 folder...\n')
        T1files=spm_select('FPList',pwd,'.*');
        [pth,nm,ext]=fileparts(T1files(1,:));                                  
        if strcmp(ext,'.nii')~=1 
        doscmd=sprintf('gdcm gdcminfo %s',deblank(T1files(1,:)));[status,cmdout]=dos(doscmd);
        %Convert enhanced dicom to standard dicom
        if strfind(cmdout,'Enhanced')
            disp('Enhanced DICOM detected. Converting to standard DICOMS...\n');
            fullfname=fullfile(t1folder,'Standard');
            if exist(fullfname,'dir')~=7
                mkdir(fullfname);
            end
            doscmd=sprintf('emf2sf --out-dir %s %s',fullfname,deblank(T1files(1,:)));
            dos(doscmd);
            t1folder=fullfname;
            cd(t1folder)
            T1files=spm_select('FPList',pwd,'.*');
        end
        hdr=spm_dicom_headers(T1files);
        spm_dicom_convertP50(hdr,'all','flat','nii');
        end
    end
    
    if ~isempty(M0folder)
        cd(M0folder); disp('In M0 folder...\n')
        M0files=spm_select('FPList',pwd,'.*');
        doscmd=sprintf('gdcm gdcminfo %s',deblank(M0files(1,:)));[status,cmdout]=dos(doscmd);
        %Convert enhanced dicom to standard dicom
        if strfind(cmdout,'Enhanced')
            display('Enhanced DICOM detected. Converting to standard DICOMS...\n');
            fullfname=fullfile(M0folder,'Standard');
            if exist(fullfname,'dir')~=7
                mkdir(fullfname);
            end
            doscmd=sprintf('emf2sf --out-dir %s %s',fullfname,deblank(M0files(1,:)));
            dos(doscmd);
            M0folder=fullfname;
            cd(M0folder)
            M0files=spm_select('FPList',pwd,'.*');
        end
        hdr=spm_dicom_headers(M0files);
        spm_dicom_convertP50(hdr,'all','flat','nii');
        
    end
    
    spm_figure('CreateWin','Interactive');
    F=spm_figure('FindWin','Interactive');set(F,'Position',[50 50 500 500]);
    spm_figure('CreateWin','Graphics');
    settings.ext='nii';
    
    
    settings.pdffolder=fullfile(settings.CBF_folder,'PDFs');
    if ~exist(settings.pdffolder,'dir')
        status=mkdir(settings.pdffolder);
    end
    
    settings.intfolder=fullfile(settings.CBF_folder,'Intermediates');
    if ~exist(settings.intfolder,'dir')
        status=mkdir(settings.intfolder);
    end
    
    [gGMCBF,SNRout,outliers,f,Slices,isquant]=asl_pipeline100317(aslfolder,t1folder,M0folder,Lesionmask,fcnl,settings,defaults);
    %     if settings.FSL==0
    %         [gGMCBF]=asl_quantification092217(Perf,M0,mask,V_CBF,settings,defaults);
    %     else
    %         %TO DO: input FSL commands
    %         %TO DO: read in FSL CBF results
    %     end
    
    %Kate: Please update location of this file.
    %fID=fopen('E:\NUdata\Swinburne/summary.txt','a+')
    fID=fopen(fullfile(aslfolder,'params.ini'),'w')
    fprintf(fID,'subtraction method=%s\n', settings.subt);
    if settings.smooth==1
        fprintf(fID,'smoothing kernel=%d\n', settings.kernel);
    else
        fprintf(fID,'smoothing kernel=0\n');
    end
    fprintf(fID,'rawPerfSNR=%4.2f\n', SNRout.rawPerfSNR);
    fprintf(fID,'temporalPerfSNR=%4.2f\n', SNRout.temporalPerfSNR);
    fprintf(fID,'rawCtrlSNR=%4.2f\n', SNRout.rawCtrlSNR);%added 070215
    fprintf(fID,'mCtrlSignal=%4.2f\n', SNRout.mCtrlSignal);%added 070215
    if length(gGMCBF)>1
        fprintf(fID,'GlobalGMCBF=%4.2f\n', gGMCBF(2));
    else
        fprintf(fID,'GlobalGMCBF=%4.2f\n', gGMCBF(1));
    end
    fprintf(fID,'deleteN=%d\n', outliers.dInd);
    fprintf(fID,'deletedPairs=%s\n', mat2str(outliers.DeletedPairs));
    fclose(fID)
    
    %% create Markdown file
    
    fileID=fopen(fullfile(aslfolder,['Report1.md']),'w');
    fprintf(fileID,'Report generated by ASLpipeline version %s in Matlab %s on %s\n',func2str(f),version, datestr(now));
    fprintf(fileID,'Please send questions & feedback to <a href=mailto:yfchen@northwestern.edu>Jennie Chen</a>.\n');
    
    fprintf(fileID,'\n111414 changed code to accomodate no T1 such that will run to completion.\n');
    fprintf(fileID,'121117 added option for LPCA DWIdenoising (Pierrick Coupe & Jose Manjon) of raw image time series.\n');
    fprintf(fileID,'070915 added option for regression of motion, 1st derivative of motion, gsignal & CSF.\n');
    fprintf(fileID,'082415 define mask using mean of ASL ctrls instead of mean, also implemented markdown for pdf generation.\n');
    fprintf(fileID,'031416 added functionality for Siemens WIP 818F BS-ASL. Modified qCBF code to remove NaN box around volume. Output indices of rejected volumes in param.ini.\n');
    fprintf(fileID,'100317 Rewrote pipeline to include denoising, PVc and normalization options.\n');
    
    fprintf(fileID,'\n##References\n');
    fprintf(fileID,'1. Alsop D., et al, MRM 2015;73:102-116\n');
    fprintf(fileID,'2. Tan H., et al, JMRI 2009;29:1134-1139\n ');
    fprintf(fileID,'3. Maumet C. et al., In: Multimodal Brain Image Analysis. MBIA 2012. Lecture Notes in Computer Science, vol 7509, pp215-224.\n ');
    fprintf(fileID,'4. Du AT., et al, Neurology 2006;67:1215-1220\n ');
    fprintf(fileID,'5. Asllani I., et al, MRM 2008;60:1362-1371\n ');
    fprintf(fileID,'6. Manjon J., et al, PLoS ONE 8(9): e73021\n ');
    fprintf(fileID,'7. Wang Z., et al, MRI 2012;30:1409-1415\n ');
    
    fprintf(fileID,'\n## ASLQA summary:\n');
    fprintf(fileID,'Name | Value\n');
    fprintf(fileID,':---------- | ----------:\n');
    if length(gGMCBF)>1
        fprintf(fileID,'Global GM CBF | %.2f ml/100g/min\n',gGMCBF(2));
    else
        fprintf(fileID,'Global GM CBF | %.2f ml/100g/min\n',gGMCBF(1));
    end
    fprintf(fileID,'rawPerfSNR | %.2f\n',SNRout.rawPerfSNR);
    fprintf(fileID,'tempPerfSNR | %.2f\n',SNRout.temporalPerfSNR);
    fprintf(fileID,'rawCtrlSNR | %.2f\n',SNRout.rawCtrlSNR);
    fprintf(fileID,'Ndeleted/Total | %d/%d\n',outliers.dInd,outliers.dInd+outliers.FinalPairNo);
    fclose(fileID);doscmd=sprintf('md2pdf Report1.md %s/Report1.pdf',settings.pdffolder);dos(doscmd)
    
    
    fileID=fopen(fullfile(aslfolder,['Report2.md']),'w');
    fprintf(fileID,'\n## ASLQA input parameters:\n');
    fprintf(fileID,'Name | Value\n');
    fprintf(fileID,':---------- | ----------:\n');
    fprintf(fileID,'FieldStrength | %d\n',settings.field);
    fprintf(fileID,'T1b | %d ms\n',settings.T1b);
    fprintf(fileID,'subtraction method | %s\n',settings.subt);
    fprintf(fileID,'Presmoothing | %d\n',settings.smooth);
    if settings.smooth~=0
        fprintf(fileID,'Presmoothing kernel | %d mm\n',settings.kernel);end
    fprintf(fileID,'fcnl | %d\n',fcnl);
    fprintf(fileID,'Outlier Removal | %s\n',settings.outliermethod);
    fprintf(fileID,'asl | %d\n',settings.asl);
    fprintf(fileID,'tau | %d\n',settings.tau);
    fprintf(fileID,'w | %d\n',settings.w);
    fprintf(fileID,'eff | %.2f\n',settings.eff);
    fprintf(fileID,'acq | %d\n',settings.acq);
    fprintf(fileID,'WIP | %d\n',settings.WIP);
    fprintf(fileID,'Regress motion&gsig | %d\n',settings.regress);
    fprintf(fileID,'denoise | %d\n',settings.NLM);
    fclose(fileID);doscmd=sprintf('md2pdf Report2.md %s/Report2.pdf',settings.pdffolder);dos(doscmd)
    
    
    fileID=fopen(fullfile(aslfolder,['Report4.md']),'w');
    fprintf(fileID,'\n## Acquisition Parameters:\n');
    fprintf(fileID,'Name | Value\n');
    fprintf(fileID,':---------- | ----------:\n');
    fprintf(fileID,'SubjectID | %s\n',settings.PatientID);
    fprintf(fileID,'SubjectAge | %s\n',settings.Age);
    fprintf(fileID,'SubjectSex | %s\n',settings.gender);
    fprintf(fileID,'Visit | %s\n',settings.visit);
    fprintf(fileID,'Series | %03d_%s\n',settings.serNumber,settings.ser);
    fprintf(fileID,'SiteID | %s\n',settings.Site);
    fprintf(fileID,'Scan Date | %s\n',settings.ScanDate);
    fprintf(fileID,'Voxel Size | %.2f %.2f %.2f\n',settings.Voxels');
    fprintf(fileID,'Slices | %d\n',Slices);
    fprintf(fileID,'TR | %d\n',settings.TR);
    fprintf(fileID,'TE | %d\n',settings.TE);
    fclose(fileID);doscmd=sprintf('md2pdf Report4.md %s/Report4.pdf',settings.pdffolder);dos(doscmd)
    
    fileID=fopen(fullfile(aslfolder,['Report5.md']),'w');
    fprintf(fileID,'![SNR](SNR.png)\n');
    fprintf(fileID,'\n## SNR summary:\n');
    fprintf(fileID,'Name | Value\n');
    fprintf(fileID,':---------- | ----------:\n');
    fprintf(fileID,'mPerfSignal | %0.2g\n',SNRout.mPerfSignal);
    fprintf(fileID,'rawPerfNoise | %0.2g\n',SNRout.rawPerfNoise);
    fprintf(fileID,'uncorr SNR | %0.2g\n',SNRout.uncorrSNR);
    fprintf(fileID,'uncorr tSNR | %0.2g\n',SNRout.uncorrtSNR);
    fprintf(fileID,'rawPerfSNR | %0.2g\n',SNRout.rawPerfSNR);
    fprintf(fileID,'temporalPerfNoise | %0.2g\n',SNRout.temporalPerfNoise);
    fprintf(fileID,'temporalPerfSNR | %0.2g\n',SNRout.temporalPerfSNR);
    fprintf(fileID,'mCtrlSignal | %0.2g\n',SNRout.mCtrlSignal);
    fprintf(fileID,'rawCtrlNoise | %0.2g\n',SNRout.rawCtrlNoise);
    fprintf(fileID,'rawCtrlSNR | %0.2g\n',SNRout.rawCtrlSNR);
    fclose(fileID);doscmd=sprintf('md2pdf Report5.md %s/Report5.pdf',settings.pdffolder);dos(doscmd)
    
    reportf=fullfile(settings.CBF_folder,sprintf('%s_%04d_%s_ASLQA.pdf',settings.visit,settings.serNumber,settings.ser));
    
    doscmd=sprintf('pdfunite %s/Report1.pdf %s/Report2.pdf %s/Report4.pdf %s/Realign.pdf',settings.pdffolder,settings.pdffolder,settings.pdffolder,settings.pdffolder);
    
    if settings.outlier==1 && strcmp(settings.outliermethod,'Zthreshold')
        doscmd=sprintf('%s %s/Report3.pdf',doscmd,settings.pdffolder);
    end
    
    doscmd=sprintf('%s %s/Report5.pdf',doscmd,settings.pdffolder);
    
    if isquant==1
        doscmd=sprintf('%s %s/globalCBF.pdf',doscmd,settings.pdffolder);
        
        if strcmp(settings.PVcmethod,'None')~=1
            doscmd=sprintf('%s %s/globalCBFPVc.pdf',doscmd,settings.pdffolder);
        end
    else
        doscmd=sprintf('%s %s/meanPerf.pdf',doscmd,settings.pdffolder);
    end
    
    doscmd=sprintf('%s %s',doscmd,reportf);
    dos(doscmd)
    
    zip(fullfile(settings.CBF_folder,'intermediates.zip'),settings.intfolder);
    rmdir(settings.intfolder,'s')
    rmdir(settings.pdffolder,'s')
    
catch err
    disp(err.message);
    disp(err.identifier);
    for k=1:length(err.stack)
        fprintf('In %s at %d\n',err.stack(k).file,err.stack(k).line);
    end
    %exit;
end
return;




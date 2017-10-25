function [gCBF,SNRout,outliers,f,Slices,isquant]=asl_pipeline100317(aslfolder,t1folder,M0folder,Lesionmask,fcnl,settings,defaults)
f=@asl_pipeline100317;
cd(aslfolder);
files=spm_select('FPList',aslfolder,['^f.*',settings.ext,'$']);
if isempty(files), files=spm_select('FPList',aslfolder,['^s.*',settings.ext,'$']);end%CYF020414

if ~isempty(M0folder)
    M0files=spm_select('FPList',M0folder,'.*nii');
    M0=M0files(1,:);clear M0files;
elseif mod(size(files,1),2)~=0%odd # of files
    M0=files(1,:);files(1,:)=[];
end

scrsz = get(0,'ScreenSize');

if ~isempty(t1folder)
    T1f=spm_vol(spm_select('FPList',t1folder,['^s.*',settings.ext]));
    T1=T1f(1,:);
    copyfile(T1.fname,fullfile(settings.CBF_folder,'T1.nii'),'f');
    copyfile(T1.fname,fullfile(settings.intfolder,'T1.nii'),'f');
    T1=fullfile(settings.intfolder,'T1.nii');
    
    disp('Coregistering first image to T1');
    fg=spm_figure('FindWin','Graphics');set(fg,'PaperPositionMode','auto','PaperType','usletter');
    x = spm_coreg(T1,files(1,:));%coregister 1st Ctrl image to T1
    %apply matrix to all images
    M  = inv(spm_matrix(x));
    MM = zeros(4,4);
    for m=1:size(files,1)
        MM = spm_get_space(files(m,:));
        spm_get_space(files(m,:), M*MM(:,:));
    end
    saveas(fg,fullfile(settings.pdffolder,'Coregister.pdf'));
    if exist('M0','var')
        x=spm_coreg(T1,M0);
        M  = inv(spm_matrix(x));
        MM = spm_get_space(M0);
        spm_get_space(M0, M*MM(:,:));
    end
    
    if ~isempty(Lesionmask)
        if ~exist(fullfile(pwd,'tmp'),'dir')
            mkdir(fullfile(pwd,'tmp'));
        end
        filenames=unzip(Lesionmask,fullfile(pwd,'tmp'));
        test=cellfun(@isempty,strfind(filenames,'gz'));
        if ~isempty(test(test==0))
            gfilenames=gunzip(filenames{test==0}); 
            
        end
        lesionT1=filenames{cellfun(@isempty,strfind(filenames,'T1.nii'))==0};
        if isempty(lesionT1), lesionT1=filenames{cellfun(@isempty,strfind(gfilenames,'T1.nii'))==0};end
        lesionmask=filenames{cellfun(@isempty,strfind(filenames,'lesionmask.nii'))==0};
        if isempty(lesionmask), lesionmask=filenames{cellfun(@isempty,strfind(gfilenames,'lesionmask.nii'))==0};end
        [pth,nm,ext]=fileparts(lesionmask);
        x=spm_coreg(T1,lesionT1);
        M=inv(spm_matrix(x));
        MM = spm_get_space(lesionmask);
        spm_get_space(lesionmask, M*MM(:,:));
        MM = spm_get_space(lesionT1);
        spm_get_space(lesionT1, M*MM(:,:));
        
        spm_reslice(char(files(1,:),lesionmask),struct('interp',1,'mask',1,'mean',0,'which',1,'wrap',[0 0 0]'));
        movefile(fullfile(pth,['r',nm,ext]),fullfile(pth,['r',nm,'_EPI',ext]));
        spm_reslice(char(T1,lesionmask),struct('interp',1,'mask',1,'mean',0,'which',1,'wrap',[0 0 0]'));
        movefile(fullfile(pth,['r',nm,ext]),fullfile(pth,['r',nm,'_T1',ext]));
        lesionmask=fullfile(pth,'rlesionmask_T1.nii');
        %remove lesion from head & segment
        headV=spm_vol(T1);
        headIM=spm_read_vols(headV);
        lesionIM=spm_read_vols(spm_vol(lesionmask));
        headIM(lesionIM>0)=0;
        headV.fname=fullfile(settings.CBF_folder,'T1_lesmasked.nii');headV.descrip='T1 with lesion masked out';
        spm_create_vol(headV);spm_write_vol(headV,headIM);
        copyfile(headV.fname,settings.intfolder);
        T1=fullfile(settings.intfolder,'T1_lesmasked.nii');
        
        lesmask=spm_read_vols(spm_vol(spm_select('FPList',pth,'^rlesionmask_EPI.nii$')));
    end
    
    %segment lesionmasked T1
    if isempty(spm_select('FPList',settings.intfolder,'^p.*T1.*nii'))
        load(which('VBMsegment.mat'));
        matlabbatch{1}.spm.tools.vbm8.estwrite.data{1} = T1;
        matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.dartelwarp.normhigh.darteltpm   = {which('Template_1_IXI550_MNI152.nii')};%070815 changed
        matlabbatch{1}.spm.tools.vbm8.estwrite.opts.tpm       = {which('TPM.nii')};%070815 changed
        try
            spm_jobman('initcfg');
            spm_jobman('run_nogui',matlabbatch);
        end
    end
    
    
    tissuefile=spm_select('FPList',settings.intfolder,'^p.*T1.*nii');
    spm_reslice(char(files(1,:),tissuefile),struct('interp',1,'mask',1,'mean',0,'which',1,'wrap',[0 0 0]'));
    tissuefile=spm_select('FPList',settings.intfolder,'^rp.*T1.*nii');
    GM=spm_read_vols(spm_vol(tissuefile(1,:)));WM=spm_read_vols(spm_vol(tissuefile(2,:)));
    CSF=spm_read_vols(spm_vol(tissuefile(3,:)));
    mask=zeros(size(GM));mask(((GM+WM+CSF)>0.1))=1;
    
end

%Motion correction - regressing out ctrl-tag pattern from motion parameters
if isempty(strfind(settings.ser,'MoCo'))%realign ONLY IF not using SMS mocoseries
    if isempty(spm_select('List',aslfolder,[settings.PatientID '_' settings.ser '_mocorrAll.png']))
        fg=spm_figure('FindWin','Graphics');set(fg,'PaperPositionMode','auto','PaperType','usletter');
        spm_realign_aslJC(files);%realign all images to 1st
        saveas(fg,fullfile(settings.pdffolder,'Realign.pdf'));
    end
end

if exist('M0','var')
    spm_reslice(char(files,M0),settings.reslcflags);
    [pth,nm,ext]=fileparts(M0);
    M0=spm_select('FPList',pth,['^r.*',settings.ext,'$']);M0=spm_read_vols(spm_vol(M0));
else
    spm_reslice(files,settings.reslcflags);
end

files=spm_select('FPList',aslfolder,['^r.*',settings.ext,'$']);
if isempty(files), files=spm_select('FPList',aslfolder,['^rs.*',settings.ext,'$']);end%CYF020414
if mod(size(files,1),2)~=0%odd # of files
    files(1,:)=[];
end

if settings.smooth==1
    disp(sprintf('Smoothing...kernel size=%d\n',settings.kernel));
    for z=1:size(files,1)
        [path name ext] = fileparts(files(z,:));
        spm_smooth(files(z,:),fullfile(path,['s' name ext]),settings.kernel);
    end
    files=spm_select('FPList',aslfolder,['^sr.*',settings.ext,'$']);
    if isempty(files), files=spm_select('FPList',aslfolder,['^f.*',settings.ext,'$']);end
    if isempty(files), files=spm_select('FPList',aslfolder,['^s.*',settings.ext,'$']);end%CYF020414
    if exist('M0','var')
        [path name ext] = fileparts(M0);
        spm_smooth(M0,fullfile(path,['s' name ext]),settings.kernel);
        M0=spm_select('FPList',M0folder,['^s.*',settings.ext,'$']);M0=spm_read_vols(spm_vol(M0(1,:)));
    end
end

if mod(size(files,1),2)~=0%odd # of files
    M0=files(1,:);files(1,:)=[];
    M0=spm_read_vols(spm_vol(M0));
end

data=spm_read_vols(spm_vol(files));ss=size(data);fac=ceil(sqrt(ss(3)));Slices=ss(3);
V_CBF=spm_vol(files(1,:));V_CBF.dt=[16 0];
V_tSNR=V_CBF;V_SNR=V_CBF;
if fcnl==1
    V_CBFser=V_CBF;V_BOLD=V_CBFser;
end
V_mask=V_CBF;
V_CBF=V_CBF;


%for determining mask
Im=mean(data(:,:,:,2:2:end),4);
thmask=zeros(size(Im));settings.thres=0.1*max(Im(:));
thmask(find(Im>settings.thres))=1;
if exist('mask','var')
    mask=mask.*thmask;
else
    mask=thmask;
end

FOVmask=mask;
if exist('lesmask','var')
    FOVmask(lesmask>0)=0;
end

V_mask.dt=[2 0];
V_mask.fname=fullfile(settings.CBF_folder,'FOVmask.nii');
spm_create_vol(V_mask);spm_write_vol(V_mask,FOVmask);

noisemask=zeros(size(mask));
noisemask([1:5 size(mask,2)-5:size(mask,2)],[1:10 size(mask,1)-10:size(mask,1)],:)=1;

paramfile=spm_select('FPList',aslfolder,['^rp_.*txt$']);
if ~isempty(paramfile)
    param=load(paramfile);param(:,1:6)=[];
    %***need to reorder JHU data or add to the following code
    if strncmpi(settings.Site,'JHU',3)==1
        temp=param;
        param(1:2:end,:)=temp(1:size(temp,1)/2,:);
        param(2:2:end,:)=temp(size(temp,1)/2+1:end,:);
    end
    param(:,4:6)=param(:,4:6)*180/pi;
    if strcmp(settings.subt,'pairwise')
            diffpara=abs(param(1:2:end,:)-param(2:2:end,:));
    else
            for a=2:size(param,1)-1
                diffpara(a-1,:)=param(a,:)-0.5*(param(a-1,:)+param(a+1,:));
            end
    end
    rawfilter=sum([abs(diffpara(:,1:3))>0.8 abs(diffpara(:,4:6))>0.8],2)>0;
    gcbf=zeros(ss(4),1);csfsig=zeros(ss(4),1);
    
    for n=1:ss(4)
        gcbf(n,1)=mean(nonzeros(data(:,:,:,n).*mask));
    end
    
    %param = [param [0 0 0 0 0 0; diff(param)]];
    
    if settings.regress==1
        param=[param gcbf];
        imtype=zeros(size(param,1),1);imtype(1:2:end)=-1;imtype(2:2:end)=1;
        covs=zeros(size(param));
        %mean centering
        covs=param-repmat(mean(param,1),[ss(4) 1]);
        %orthogonalizing
        for n=1:size(param,2)
            covs(:,n)=covs(:,n)-imtype*regress(covs(:,n),imtype);end
        
        rdata=reshape(data,prod(ss(1:3)),ss(4))';
        ndata =rdata - (covs * inv(covs' * covs)* covs') * rdata;
        ndata =reshape(ndata', ss);
        clear rdata covs;
    else
        ndata=data;
    end
else
    ndata=data;
    rawfilter=zeros(length(gcbf),1);
end

Ctrl=ndata(:,:,:,2:2:end);
Tag=ndata(:,:,:,1:2:end);
if isempty(M0folder)
    if exist('M0','var')==0 && settings.is3D==0
        M0=Ctrl;
    end
end


if size(Ctrl,4)>size(Tag,4)
    msg=sprintf('Warning: Ctrl volumes=%g > Tag volumes=%g!',size(Ctrl,4),size(Tag,4));
    display(msg);
    Ctrl(:,:,:,size(Ctrl,4)-1)=[];
elseif size(Tag,4)>size(Ctrl,4)
    msg=sprintf('Warning: Tag volumes=%g > Ctrl volumes=%g!',size(Tag,4),size(Ctrl,4));
    display(msg);
    Tag(:,:,:,size(Tag,4)-1)=[];
end

V_CBF.fname=fullfile(settings.CBF_folder,'meanEPI.nii');
spm_create_vol(V_CBF);spm_write_vol(V_CBF,mean(Ctrl,4).*mask);

if strcmp(settings.subt,'pairwise')
        Perf=Ctrl-Tag;SubOrder='Ctrl-Tag';
        if fcnl==1, BOLD=(Ctrl+Tag)/2;end
        subtmethod='pair';
else
        SubOrder='Ctrl-Tag';
        for a=1:size(Tag,4)-1
            Perf(:,:,:,a*2-1)=Ctrl(:,:,:,a)-(Tag(:,:,:,a)+Tag(:,:,:,a+1))/2;
            Perf(:,:,:,a*2)=(Ctrl(:,:,:,a)+Ctrl(:,:,:,a+1))/2-Tag(:,:,:,a+1);
            if settings.asl==2&&isempty(M0folder)
                M0(:,:,:,a*2-1)=Ctrl(:,:,:,a);
                M0(:,:,:,a*2)=(Ctrl(:,:,:,a)+Ctrl(:,:,:,a+1))/2;
            end
            if fcnl==1
                BOLD(:,:,:,a*2-1)=(Ctrl(:,:,:,a)+(Tag(:,:,:,a)+Tag(:,:,:,a+1))/2)/2;
                BOLD(:,:,:,a*2)=((Ctrl(:,:,:,a)+Ctrl(:,:,:,a+1))/2+Tag(:,:,:,a))/2;
            end
        end
        subtmethod='surr';%JC added 3/8/13
end

if mean(nonzeros(mean(Perf,4).*mask))<0
    Perf=-Perf;
    SubOrder='Tag-Ctrl';
end
Perf(isnan(Perf))=0;

%calculate uncorrected SNR and tSNR
meanuncorrPerf=mean(Perf,4);
SNRmap=meanuncorrPerf;
for slc=1:Slices
    acqcorr=(slc-1)*settings.acq;
    bkg(slc,1)=std(nonzeros(noisemask(:,:,slc).*meanuncorrPerf(:,:,slc)));
    SNRmap(:,:,slc)=SNRmap(:,:,slc).*mask(:,:,slc)/(exp(-acqcorr/settings.T1b)*bkg(slc,1));
end

SNRmap(isnan(SNRmap)|isinf(SNRmap))=0;
uncorrSNR=mean(nonzeros(SNRmap));clear SNRmap meanuncorrPerf bkg;

noise=std(Perf.*repmat(mask,[1 1 1 size(Perf,4)]),1,4);noise(isnan(noise))=0;
tSNR=mask.*(mean(Perf.*repmat(mask,[1 1 1 size(Perf,4)]),4)./noise);
tSNR(isnan(tSNR))=0;tSNR(isinf(tSNR))=0;
uncorrtSNR=mean(nonzeros(tSNR));clear noiss tSNR;


outliers=struct('dInd',[],'DeletedPairs',0,'FinalPairNo',size(Perf,4)); 
delInd=zeros(size(Perf,4),1);

if settings.outlier==1
    switch settings.outliermethod
        case 'Robust'
            disp('Robust outlier removal using M-estimator.');
            %output Perf maps for M-Estimator
            for a=1:size(Perf,4)
                %correct for slice acq effects
                for slc=1:Slices
                    acqcorr=(slc-1)*settings.acq;
                    Perf(:,:,slc,a)=Perf(:,:,slc,a)/(exp(-acqcorr/settings.T1b));end
                V_CBF.fname=fullfile(settings.CBF_folder,sprintf('tmpPerf_%04d.%s',a,settings.ext));
                spm_create_vol(V_CBF);spm_write_vol(V_CBF,Perf(:,:,:,a));
            end
            inputFiles=spm_select('FPList',settings.CBF_folder,['^tmp.*',settings.ext,'$']);
            asl_robustJC(cellstr(inputFiles), fullfile(settings.intfolder,sprintf('%s_RobustPerf.%s',settings.PatientID,settings.ext)), 1.345);
            delete(fullfile(settings.CBF_folder,'tmp*'));
            
            robustfile=spm_select('FPList',settings.intfolder,sprintf('%s_RobustPerf.%s',settings.PatientID,settings.ext));
            outliers.dInd=length(find(delInd));
            outliers.DeletedPairs=find(delInd);
            outliers.FinalPairNo=size(Perf,4)-length(find(delInd));
            Perf=spm_read_vols(spm_vol(robustfile));
            
            asl_robustJC(cellstr(files(2:2:end,:)), fullfile(settings.intfolder,sprintf('%s_RobustM0.%s',settings.PatientID,settings.ext)), 1.345);
            M0=spm_read_vols(spm_vol(fullfile(settings.intfolder,sprintf('%s_RobustM0.%s',settings.PatientID,settings.ext))));
        case 'Zthreshold'
            disp('Outlier removal using Z-threshold.');
            for a=1:size(Perf,4)
                %correct for slice acq effects
                for slc=1:Slices
                    acqcorr=(slc-1)*settings.acq;
                    Perf(:,:,slc,a)=Perf(:,:,slc,a)/(exp(-acqcorr/settings.T1b));end
                meanval(a,1)=mean(nonzeros(mask.*Perf(:,:,:,a)));
                stdevval(a,1)=std(nonzeros(mask.*Perf(:,:,:,a)));
            end
            meanfilt=zeros(size(meanval));stdevfilt=meanfilt;
            if log10(max(stdevval)-min(stdevval))>1
                meanfilt(((meanval>(mean(meanval)+2.5*std(meanval)))|(meanval<(mean(meanval)-2.5*std(meanval)))))=1;
                stdevfilt((stdevval>(mean(stdevval)+1.5*std(stdevval))))=1;
            end
            
            fg=figure('Position',[100 0 scrsz(3)*.75 scrsz(4)/2],'Color','white','PaperPositionMode','auto','PaperType','usletter');
            subplot(1,2,1);plot(meanval,'-o');hold on;plot(1:length(meanval),repmat(mean(meanval),[1 length(meanval)]),'r');
            plot(1:length(meanval),repmat(mean(meanval)+2.5*std(meanval),[1 length(meanval)]),'r--');
            plot(1:length(meanval),repmat(mean(meanval)-2.5*std(meanval),[1 length(meanval)]),'r--');title('mean brain signal');hold off
            
            subplot(1,2,2);plot(stdevval,'-o');hold on;plot(1:length(stdevval),repmat(mean(stdevval),[1 length(stdevval)]),'r');
            plot(1:length(stdevval),repmat(mean(stdevval)+1.5*std(stdevval),[1 length(stdevval)]),'r--');title('stdev brain signal');hold off
            
            subplot(1,2,1);hold on;plot(find(meanfilt==1),meanval(find(meanfilt==1)),'rx');
            subplot(1,2,2);hold on;plot(find(stdevfilt==1),stdevval(find(stdevfilt==1)),'rx');
            saveas(fg,'OutlierDetection.png');
            
            
            delInd=(rawfilter+meanfilt+stdevfilt)>0;
            PerfO=Perf;mean_PerfO=mean(PerfO,4);
            Perf(:,:,:,find(delInd))=[];
            if exist('M0','var') && size(M0,4)>1
                M0O=M0;
                M0(:,:,:,find(delInd))=[];
            end
            
            if fcnl==1
                BOLDO=BOLD;
                BOLD(:,:,:,find(delInd))=[];
            end
            outliers.dInd=length(find(delInd));
            outliers.DeletedPairs=find(delInd);
            outliers.FinalPairNo=size(PerfO,4)-length(find(delInd));
            
            fileID=fopen(fullfile(aslfolder,['Report3.md']),'w');
            fprintf(fileID,'![Outlier](OutlierDetection.png)\n');
            fprintf(fileID,'Deleted images [%s] based on |mean|>2.5*std(mean)\n',num2str(find(meanfilt==1)'));
            fprintf(fileID,'Deleted images [%s] based on stdev>1.5*std(stdev)\n',num2str(find(stdevfilt==1)'));
            fprintf(fileID,'Deleted images [%s] based on motion\n',num2str(find(rawfilter)'));
            fprintf(fileID,'Total: deleted %s/%s images: %s .\n',num2str(outliers.dInd),num2str(size(PerfO,4)),num2str(find(delInd)'));
            fclose(fileID);
            doscmd=sprintf('md2pdf %s/Report3.md %s/Report3.pdf',aslfolder,settings.pdffolder);dos(doscmd)
            
            
        case 'None'
            disp('No outlier removal.');
            outliers.dInd=length(find(delInd));
            outliers.DeletedPairs=find(delInd);
            outliers.FinalPairNo=size(Perf,4)-length(find(delInd));
            
    end
end

%added 070215--------------------------
if settings.NLM==1
    % parameters for block-wise NLM
    M=3;    %  controlling size of search volume
    alpha=1;  % controlling block size
    
    denoisedPerf=zeros(size(Perf));I_orig=Perf;
    disp('Denoising perfusion timeseries...');
    for i=1:size(Perf,4)
        ref=I_orig(:,:,:,i);
        %estimation of noise level
        estimated_noise=RicianSTD(I_orig(:,:,:,i));
        
        s=size(I_orig(:,:,:,i));
        
        % adjust the noise level to balance the effect of denoising and smoothing
        level=estimated_noise/4;
        
        % Make sure the intensity is non-negative
        L = min(ref(:));
        I_orig(:,:,:,i) = I_orig(:,:,:,i) + abs(L);
        
        % Obtaining images with preserved featrues
        I_SP(:,:,:,i)=onlm(I_orig(:,:,:,i),M,alpha,level); % onlm - optimized block-wise non-local means denoising, compiled mex-file
        I_SP(:,:,:,i)=I_SP(:,:,:,i) - abs(L);
        
        
        % Obtaining images with noise components removed
        I_BP(:,:,:,i)=onlm(I_orig(:,:,:,i),M,alpha+1,level);
        I_BP(:,:,:,i)=I_BP(:,:,:,i) - abs(L);
        
        % Further denoising with DT-CWT
        denoisedPerf(:,:,:,i) = DTCWT(I_orig(:,:,:,i),I_SP(:,:,:,i),I_BP(:,:,:,i));
    end
    disp('Done!');
    Perf=denoisedPerf;
    clear denoisedPerf I_orig I_BP I_SP L M alpha level;
end
%added 070215--------------------------

%output Perf & M0
save_avw(Perf,fullfile(settings.intfolder,sprintf('%s_%s_Perf.%s',settings.visit,settings.ser,settings.ext)),'f',[settings.Voxels' settings.TR/1000]);
if exist('M0','var')
    save_avw(M0,fullfile(settings.intfolder,sprintf('%s_%s_M0.%s',settings.visit,settings.ser,settings.ext)),'f',[settings.Voxels' settings.TR/1000]);end

if fcnl==1
    save_avw(BOLD,fullfile(settings.intfolder,sprintf('%s_%s_BOLD.%s',settings.visit,settings.ser,settings.ext)),'f',[settings.Voxels' settings.TR/1000]);
end
%SNR calculations
mean_Perf=mean(Perf,4);
SNRmap=mean_Perf;
cSNRmap=mean(Ctrl,4);
for slc=1:Slices
    y=ceil(slc/fac);
    x=slc-(y-1)*fac;
    bkg(slc,:)=[std(nonzeros(noisemask(:,:,slc).*mean_Perf(:,:,slc))) std(nonzeros(noisemask(:,:,slc).*mean(Ctrl(:,:,slc,:),4)))];
    SNRmap(:,:,slc)=(SNRmap(:,:,slc).*mask(:,:,slc))/bkg(slc,1);
    cSNRmap(:,:,slc)=(cSNRmap(:,:,slc).*mask(:,:,slc))/bkg(slc,2);
    mosSNR((y-1)*size(Ctrl,2)+1:y*size(Ctrl,2),(x-1)*size(Ctrl,1)+1:x*size(Ctrl,1))=imrotate(SNRmap(:,:,slc),90);
end

SNRmap(isnan(SNRmap)|isinf(SNRmap))=0;cSNRmap(isnan(cSNRmap)|isinf(cSNRmap))=0;
[row,col]=find(isnan(bkg));bkg(row,:)=[];
meansig=mean(nonzeros(mean_Perf.*mask));
meanCtrl=mean(nonzeros(Ctrl.*repmat(mask,[1,1,1,size(Ctrl,4)])));
fg=figure('Position',[100 0 scrsz(3)*.75 scrsz(4)/2],'Color','white','PaperPositionMode','auto','PaperType','usletter');
if strcmp(settings.outliermethod,'Robust')==1%Robust outlier detection only outputs mean_perf map
    SNRout=struct('mPerfSignal',mean(mean_Perf((mask>0))),'rawPerfNoise',mean(bkg(:,1)),'uncorrSNR',uncorrSNR,'uncorrtSNR',uncorrtSNR,...
        'rawPerfSNR',mean(SNRmap((mask>0))),'temporalPerfNoise',0,'temporalPerfSNR',NaN,...
        'mCtrlSignal',meanCtrl,'rawCtrlNoise',mean(bkg(:,2)),'rawCtrlSNR',mean(cSNRmap((mask>0))));
    imagesc(mosSNR,[-10 3*mean(nonzeros(SNRmap))]);colorbar('vert');title('Raw SNR');colormap gray;
else
    noise=std(Perf.*repmat(mask,[1 1 1 size(Perf,4)]),1,4);noise(isnan(noise))=0;
    tSNR=mask.*(mean(Perf.*repmat(mask,[1 1 1 size(Perf,4)]),4)./noise);
    tSNR(isnan(tSNR))=0;tSNR(isinf(tSNR))=0;
    
    SNRout=struct('mPerfSignal',mean(mean_Perf((mask>0))),'rawPerfNoise',mean(bkg(:,1)),'uncorrSNR',uncorrSNR,'uncorrtSNR',uncorrtSNR,...
        'rawPerfSNR',mean(SNRmap((mask>0))),'temporalPerfNoise',mean(noise((mask>0))),'temporalPerfSNR',mean(tSNR((mask>0))),...
        'mCtrlSignal',meanCtrl,'rawCtrlNoise',mean(bkg(:,2)),'rawCtrlSNR',mean(cSNRmap((mask>0))));
    for n=1:Slices
        y=ceil(n/fac);
        x=n-(y-1)*fac;
        mostSNR((y-1)*size(Ctrl,2)+1:y*size(Ctrl,2),(x-1)*size(Ctrl,1)+1:x*size(Ctrl,1))=imrotate(tSNR(:,:,n),90);
    end
    subplot(1,2,1);imagesc(mosSNR,[0 3*mean(nonzeros(SNRmap))]);colorbar('vert');title('Raw SNR');colormap gray;
    subplot(1,2,2);imagesc(mostSNR,[0 3*mean(nonzeros(tSNR))]);colorbar('vert');title('Temporal SNR');colormap gray;
    V_tSNR.fname=fullfile(settings.CBF_folder,sprintf('%s_%s_tSNR_%s.%s',settings.visit,settings.ser,subtmethod,settings.ext));
    V_tSNR.descrip=sprintf('%s_%g pair tSNR*10 image',SubOrder,size(Perf,4));
    V_tSNR=spm_create_vol(V_tSNR);
    V_tSNR=spm_write_vol(V_tSNR,tSNR*10);
end
saveas(fg,'SNR.png');

V_SNR.fname=fullfile(settings.CBF_folder,sprintf('%s_%s_rawSNR_%s.%s',settings.visit,settings.ser,subtmethod,settings.ext));
V_SNR.descrip=sprintf('rSNR map of %s_%g pair averaged perfusion image',SubOrder,size(Perf,4));
V_SNR=spm_create_vol(V_SNR);
V_SNR=spm_write_vol(V_SNR,SNRmap);

if exist('M0','var')
    gCBF=asl_quantification(Perf,M0,mask,V_CBF,settings,defaults);
    isquant=1;
else
    gCBF=mean(mean_Perf(mask>0));
    isquant=0;
    disp('WARNING: no M0, so no quantification!')
    for n=1:size(mean_Perf,3)
    y=ceil(n/fac);
    x=n-(y-1)*fac;
    temp1((y-1)*size(mean_Perf,2)+1:y*size(mean_Perf,2),(x-1)*size(mean_Perf,1)+1:x*size(mean_Perf,1))=imrotate(mean_Perf(:,:,n).*mask(:,:,n),90);
end
fg=figure('Position',[500 scrsz(4) scrsz(3)*.25 scrsz(4)/2],'color','white','PaperPositionMode','auto','PaperType','usletter');
subplot(2,1,1);imagesc(temp1,[0 120]);colormap(jet);colorbar('vert');
title('mean Perfusion map')
subplot(2,1,2);hist(mean_Perf(mask>0),128);xlabel('mean Perf signal (a.u.)');ylabel('# voxels');title('Histogram of mean perfusion signals');
saveas(fg,fullfile(settings.pdffolder,'meanPerf.pdf'));
end

mat=spm_select('FPList',settings.intfolder,'^y.*nii$');
T1=spm_select('FPList',settings.intfolder,'^T1.*nii$');
tissuefile=spm_select('FPList',settings.intfolder,'^p.*nii$');
mask=spm_select('FPList',settings.CBF_folder,'^FOVmask.nii$');
CBFfiles=spm_select('FPList',settings.CBF_folder,'.*qCBF.*nii$');
meanEPI=spm_select('FPList',settings.CBF_folder,'^meanEPI.*nii$');
if ~isempty(mat)
    load(which('vbmwarp.mat'));
    IMfiles=char(CBFfiles,mask,meanEPI);
    if exist('T1','var')
        IMfiles=char(IMfiles,T1,tissuefile);
    end
    matlabbatch{1,1}.spm.tools.vbm8.tools.defs.images=cellstr(IMfiles);
    matlabbatch{1,1}.spm.tools.vbm8.tools.defs.field1=cellstr(mat);
    
    try
        spm_jobman('initcfg');
        spm_jobman('run_nogui',matlabbatch);
    end
end
end

%Modification Notes
%092217 Combined quantification and two PVc techniques
%092217 To Do: incorporate RobustPerf quantification/PVc
function [gCBF,settings]=asl_quantification(Perf,M0,mask,V_qCBF,settings,defaults)
f=@asl_quantification;
scrsz = get(0,'ScreenSize');

%Quantification
mean_Perf=mean(Perf,4);
buf=zeros(size(mean_Perf));
Slice=size(Perf,3);
fac=ceil(sqrt(Slice));
if settings.asl==1%PASL
    buf((mask>0))=(100*0.9*mean_Perf((mask>0)))./(2*settings.eff*M0((mask>0))*(settings.tau/60000)*exp(-settings.w/settings.T1b));
else%CASL/pCASL
    buf((mask>0))=(100*0.9*mean_Perf((mask>0))./(2*settings.eff*M0((mask>0))*(settings.T1b/60000)*(exp(-settings.w/settings.T1b)-exp(-(settings.tau+settings.w)/settings.T1b))));
end


if settings.WIP==1
    buf=buf/10;
end

qCBF=mean(buf,4);
qCBF((qCBF<-10|qCBF>150))=0;buf((buf<-10|buf>150))=0;
V_qCBF.fname=fullfile(settings.CBF_folder,sprintf('%s_%04d_qCBF.%s',settings.visit,settings.serNumber,settings.ext));
V_qCBF.descrip=sprintf('qCBFmap');
V_qCBF=spm_create_vol(V_qCBF);
V_qCBF=spm_write_vol(V_qCBF,qCBF);

gCBF=mean(nonzeros(qCBF));

%Code for displaying qCBFmaps--------------------------
for n=1:size(qCBF,3)
    y=ceil(n/fac);
    x=n-(y-1)*fac;
    temp1((y-1)*size(qCBF,2)+1:y*size(qCBF,2),(x-1)*size(qCBF,1)+1:x*size(qCBF,1))=imrotate(qCBF(:,:,n).*mask(:,:,n),90);
end
fg=figure('Position',[500 scrsz(4) scrsz(3)*.25 scrsz(4)/2],'color','white','PaperPositionMode','auto','PaperType','usletter');
subplot(2,1,1);imagesc(temp1,[0 120]);colormap(jet);colorbar('vert');
title(['global CBF\_uncorr=',sprintf('%.2f ml/100g/min',gCBF)])
subplot(2,1,2);hist(nonzeros(qCBF),128);xlabel('qCBF (ml/100g/min)');ylabel('# voxels');title('Histogram of uncorr perfusion values');
saveas(fg,fullfile(settings.pdffolder,'globalCBF.pdf'));

if strcmp(settings.PVcmethod,'none')~=1 
    gCBF=PVcorrection(Perf, M0, settings);
end
end

function gGMCBF=PVcorrection(Perf,M0,settings)

%locate the relevant files
qCBFf=spm_select('FPList',settings.CBF_folder,'.*qCBF.nii$');qCBFf=qCBFf(1,:);
V_qCBF=spm_vol(qCBFf);
qCBF=spm_read_vols(V_qCBF);

meanEPI=spm_select('FPList',settings.CBF_folder,'^meanEPI.nii$');
mask=spm_read_vols(spm_vol(spm_select('FPList',settings.CBF_folder,'^FOVmask.nii$')));
head=spm_select('FPList',settings.CBF_folder,'^T1.nii$');

tissuefile=spm_select('FPList',settings.intfolder,'^rp.*T1.*nii');
GM=spm_read_vols(spm_vol(tissuefile(1,:)));WM=spm_read_vols(spm_vol(tissuefile(2,:)));
CSF=spm_read_vols(spm_vol(tissuefile(3,:)));
GM(mask<1)=0;WM(mask<1)=0;CSF(mask<1)=0;

switch settings.PVcmethod
    case 'Asllani'
        [GMdm, WMdm]=pvc_Asllani(Perf,GM,WM,CSF,mask,5,0);
        [GMM0, WMM0]=pvc_Asllani(M0,GM,WM,CSF,mask,5,1);
        
        qCBF_PVc=zeros(size(GMdm));tmpGM=qCBF_PVc;tmpWM=qCBF_PVc;
        if settings.asl==1%PASL
            tmpGM=100*0.98*mask.*mean(GMdm,4)./(2*settings.eff*mean(GMM0,4)*(settings.tau/60000)*exp(-settings.w/settings.T1b));
            tmpWM=100*0.82*mask.*mean(WMdm,4)./(2*settings.eff*mean(WMM0,4)*(settings.tau/60000)*exp(-settings.w/settings.T1b));
        else%CASL/pCASL
            tmpGM=100*0.98*mask.*mean(GMdm,4)./(2*settings.eff*mean(GMM0,4)*(settings.T1b/60000)*(exp(-settings.w/settings.T1b)-exp(-(settings.tau+settings.w)/settings.T1b)));
            tmpWM=100*0.82*mask.*mean(WMdm,4)./(2*settings.eff*mean(WMM0,4)*(settings.T1b/60000)*(exp(-settings.w/settings.T1b)-exp(-(settings.tau+settings.w)/settings.T1b)));
        end
        
        tmpGM(GM<0.01)=0;tmpWM(WM<0.01)=0;
        if settings.WIP==1, tmpGM=tmpGM/10;tmpWM=tmpWM/10;end
        V_qCBF.fname=fullfile(settings.CBF_folder,sprintf('%s_%04d_GMqCBF_LRPVc.%s',settings.visit,settings.serNumber,settings.ext));
        V_qCBF.descrip=sprintf('Pure GM qCBFmap with LR(Asllani) PVc');
        V_qCBF=spm_create_vol(V_qCBF);
        V_qCBF=spm_write_vol(V_qCBF,tmpGM);
        
        V_qCBF.fname=fullfile(settings.CBF_folder,sprintf('%s_%04d_WMqCBF_LRPVc.%s',settings.visit,settings.serNumber,settings.ext));
        V_qCBF.descrip=sprintf('Pure WM qCBFmap with LR(Asllani) PVc');
        V_qCBF=spm_create_vol(V_qCBF);
        V_qCBF=spm_write_vol(V_qCBF,tmpWM);
        
        gCBF(2)=mean(nonzeros(tmpGM));
        qCBF_PVc=tmpGM;
        CBFfiles=spm_select('FPlist',settings.CBF_folder,'.*LRPVc.nii');
    case 'PET'
        GM((GM<0.3))=0;
        qCBF_PVc=zeros(size(qCBF));
        qCBF_PVc((GM>0))=qCBF((GM>0))./(GM((GM>0))+0.4*WM((GM>0)));
        qCBF_PVc(qCBF_PVc<0|qCBF_PVc>200)=0;
        gCBF(2)=mean(qCBF_PVc((GM>0)));
        V_qCBF.fname=fullfile(settings.CBF_folder,sprintf('%s_%04d_qCBF_PETPVc.%s',settings.visit,settings.serNumber,settings.ext));
        V_qCBF.descrip=sprintf('qCBFmap with PET PVc');
        V_qCBF=spm_create_vol(V_qCBF);
        V_qCBF=spm_write_vol(V_qCBF,qCBF_PVc);
        CBFfiles=spm_select('FPlist',CBF_folder,'.*PETPVc.nii');
end

fac=ceil(sqrt(size(qCBF,3)));
for n=1:size(qCBF,3)
    y=ceil(n/fac);
    x=n-(y-1)*fac;
    temp1((y-1)*size(qCBF,2)+1:y*size(qCBF,2),(x-1)*size(qCBF,1)+1:x*size(qCBF,1))=imrotate(qCBF(:,:,n).*mask(:,:,n),90);
    temp2((y-1)*size(qCBF_PVc,2)+1:y*size(qCBF_PVc,2),(x-1)*size(qCBF_PVc,1)+1:x*size(qCBF_PVc,1))=imrotate(qCBF_PVc(:,:,n).*mask(:,:,n),90);
end
V=spm_vol(qCBFf);
spm_create_vol(V);spm_write_vol(V,qCBF.*mask);

gGMCBF(1)=mean(qCBF(GM>0));
gGMCBF(2)=mean(qCBF_PVc((GM>0)));
scrsz = get(0,'ScreenSize');
fg=figure('Position',[500 scrsz(4)/2 scrsz(3)*.5 scrsz(4)-100],'color','white','PaperPositionMode','auto','PaperType','usletter');
subplot(2,2,1);imagesc(temp1,[0 150]);colormap(jet);colorbar('vert'); title(['GMCBF\_uncorr=',sprintf('%.2f ml/100g/min',gGMCBF(1))])
subplot(2,2,2);imagesc(temp2,[0 150]);colormap(jet);colorbar('vert');title(['GMCBF\_PVc=',sprintf('%.2f ml/100g/min',gGMCBF(2))])
subplot(2,2,3);hist(nonzeros(qCBF),128);xlabel('qCBF (ml/100g/min)');ylabel('# voxels');title('Histogram of uncorr perfusion values');
subplot(2,2,4);hist(nonzeros(qCBF_PVc),128);xlabel('qCBF (ml/100g/min)');ylabel('# voxels');title('Histogram of PVc perfusion values');
saveas(fg,fullfile(settings.pdffolder,'globalCBFPVc.pdf'));

end


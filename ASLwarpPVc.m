function ASLwarpPVc(smri_directory, asl_base)
disp(sprintf('Running ASL_RobustNormalization pipeline, date=%s...\n',datestr(now)))
disp(sprintf('SMRIdir=%s...\n',smri_directory))
disp(sprintf('ASLdir=%s...\n',asl_base))
try

for resource= dir(sprintf('%s/sequence*',asl_base))
  ASLoutput = sprintf('%s/%s',asl_base,resource.name);


%locate the files in ASLoutput
CBF=spm_select('FPList',deblank(ASLoutput),'.*qCBF.*nii$');CBF=CBF(1,:);
T1ASL=spm_select('FPList',deblank(ASLoutput),'^T1.nii$');
%locate desired files in RobustOUtput
field=spm_select('FPList',smri_directory,'anat2tpl.warp.field.nii$');%replace smri_directory with temporary folder in which the tar file contents are located
head=spm_select('FPList',smri_directory,'^head.nii$');%replace smri_directory with temporary folder in which the tar file contents are located
others=cellstr(strvcat(spm_select('FPList',smri_directory,'^brain.msk.nii$'),...
    spm_select('FPList',smri_directory,'^gm.nii$'),...
    spm_select('FPList',smri_directory,'^wm.nii$')));%replace smri_directory with temporary folder in which the tar file contents are located
%load matlabbatch file and perform 
%1) coregistration of T1ASL to head.nii
%2) warping all files to template space using a nat2tpl.warp.field.nii
%load('C:\Users\Admin\Documents\MATLAB\NUNDA\warp.mat');%replace with NUNDA location
load('/projects/p20394/software/matlab/coreg.mat');
matlabbatch{1,1}.spm.spatial.coreg.estimate.ref{1}=head;
matlabbatch{1,1}.spm.spatial.coreg.estimate.source{1}=T1ASL;
matlabbatch{1,1}.spm.spatial.coreg.estimate.other=cellstr(CBF);
test=load('/projects/p20394/software/matlab/warp.mat');
matlabbatch{1,2}=test.matlabbatch{1,1};
matlabbatch{1,2}.spm.tools.vbm8.tools.defs.field1{1}=field;
matlabbatch{1,2}.spm.tools.vbm8.tools.defs.images=cat(1,head,others,cellstr(CBF));
try
    spm_jobman('initcfg');
    spm_jobman('run_nogui',matlabbatch);
end

%Locate the warped template space tissue & brain masks
GM=spm_read_vols(spm_vol(spm_select('FPList',smri_directory,'^wgm.nii$')));%replace smri_directory with temporary folder in which the tar file contents are located
WM=spm_read_vols(spm_vol(spm_select('FPList',smri_directory,'^wwm.nii$')));%replace smri_directory with temporary folder in which the tar file contents are located
brain=spm_read_vols(spm_vol(spm_select('FPList',smri_directory,'^wbrain.msk.nii$')));%replace smri_directory with temporary folder in which the tar file contents are located
brain(isnan(brain))=0;
%located the template space ASL file and mask to >5
ASLf=cellstr(spm_select('FPList',deblank(ASLoutput),'^w.*qCBF.*nii$'));
test=cellfun(@isempty,strfind(ASLf,'_PVc'));ASLf(find(test==0))=[];
ASL=spm_read_vols(spm_vol(char(ASLf)));
ASL(isnan(ASL))=0;
ASL(find(ASL<5))=0;
%Smooth tissue files to native resolution of ASL to match degree of
%smoothness
sGM=zeros(size(GM));spm_smooth(GM,sGM, [3.44 3.44 4]);
sWM=zeros(size(WM));spm_smooth(WM,sWM, [3.44 3.44 4]);
sWM99=sWM;sWM99(find(sWM99<0.99))=0;sWM99(find(sWM99))=1;
%Limit calculations to GM>=0.3 to minimize elevated values due to division
%by small numbers
sGM(find(sGM<0.3))=0;
V=spm_vol(spm_select('FPList',smri_directory,'wgm.nii$'));
[pth,nm,ext]=fileparts(V.fname);
V.fname=fullfile(char(ASLoutput),['wsGM',ext]);
V.descrip='smoothed GM probability > 0.3';
spm_create_vol(V);spm_write_vol(V,sGM);
V.fname=fullfile(char(ASLoutput),['wsWM99',ext]);
V.descrip='smoothed WM mask > 0.99';
spm_create_vol(V);spm_write_vol(V,sWM99);

%Define FOVmask, which combines (brain-lesion), GM>0.3 and ASL>5
FOV=zeros(size(ASL));GMmask=sGM;GMmask(find(sGM))=1;
FOV=GMmask.*brain;FOV(find(ASL<5))=0;
V=spm_vol(spm_select('FPList',smri_directory,'^wbrain.msk.nii$'));
[pth,nm,ext]=fileparts(V.fname);
V.fname=fullfile(char(ASLoutput),['FOVmask',ext]);spm_create_vol(V);spm_write_vol(V,FOV);
%Perform partial volume correction by subtracting Pwm*WM_CBF and then
%division by Pgm
WMCBF=mean(ASL(find(sWM99)));
ASL_PV=ASL-sWM*WMCBF;
ASL_PV(find(sGM<0.3))=0;ASL_PV(find(sGM>=0.3))=ASL_PV(find(sGM>=.3))./sGM(find(sGM>=0.3));
V=spm_vol(char(ASLf));
[pth,nm,ext]=fileparts(V.fname);
V.fname=fullfile(pth,[nm,'_PVc',ext]);V.descrip='Partial volume corrected (GM>=0.3) qCBF ';
spm_create_vol(V);spm_write_vol(V,ASL_PV);

% ASL_PV=ASL;ASL_PV(find(sGM<0.3))=0;
% ASL_PV(find(sGM>=0.3))=ASL_PV(find(sGM>=.3))./(sGM(find(sGM>=0.3))+0.4*sWM(find(sGM>=0.3)));
% V.fname=fullfile(pth,[nm,'_PVcPET',ext]);V.descrip='Partial volume corrected (GM>=0.3) qCBF based on GM:WM=2.5';

end

catch err
    disp(err.message);
    disp(err.identifier);
    for k=1:length(err.stack)
        fprintf('In %s at %d\n',err.stack(k).file,err.stack(k).line);
    end
    %exit;
end
return;

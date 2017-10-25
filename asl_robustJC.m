% Compute a robust perfusion-weigthed map using Huber's M-estimator from a
% series of perfusion-weithed volumes (repeated pair-subtracted ASL
% volumes)
%
% Dependencies: SPM8 or more recent
%
% inputFiles: Path to the repeated pair-substracted ASL volumes as a
% cell-array
% outputFile: Path to the output file
% huberParam: Parameter of Huber's M-estimator (1.5 by default)
%
% if inputFiles or outputFile is omitted then the parameters are asked
% interactively using SPM batch GUI. If both are defined the computation is
% performed without user interaction.
%_______________________________________________________________________
% Copyright (C) 2013 Camille Maumet
function asl_robustJC(inputFiles, outputFile, huberParam)
    % By default set Huber's M-estimator parameter to 1.5
    if ~exist('huberParam', 'var')
        huberParam = 1.5;
    end

    % Define an imcalc job
    matlabbatch{1}.spm.util.imcalc.expression = ...
     ['cellfun(@(y) robust_exclude_nan(y, ' num2str(huberParam) '), num2cell(X,1));'];
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
    
    if nargin > 0
        inputImg = nifti(inputFiles);
        % 4-D nifti volume
        if numel(inputImg.dat.dim) == 4
            numVolumes = inputImg.dat.dim(4);
            for i = 1:numVolumes
                matlabbatch{1}.spm.util.imcalc.input{i,1} = [inputFiles ',' num2str(i)];
            end
        else
            matlabbatch{1}.spm.util.imcalc.input = inputFiles;
        end
    end
    if nargin > 1
        [outputDir, outputName, outputExt] = spm_fileparts(outputFile);
        matlabbatch{1}.spm.util.imcalc.output = [outputName outputExt];
        matlabbatch{1}.spm.util.imcalc.outdir = {outputDir};
        
        matlabbatch{1}.spm.util.imcalc.output = outputFile;
        spm_jobman('run', matlabbatch);
    else
        spm_jobman('interactive', matlabbatch);
    end
end
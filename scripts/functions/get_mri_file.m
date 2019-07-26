function [mri_file, update_resl, resl_new] = get_mri_file(subj_id, mri_dir)
    % Get the right file when loading MRIs, stores long filenames.
    switch subj_id
        case 's0001'
            mri_file      = fullfile(mri_dir, 'filename.nii');
        case '0002'
            mri_file      = fullfile(mri_dir, 'filename.nii');
        case '0003'
            mri_file      = fullfile(mri_dir, 'filename.nii');
        case '0004'
            mri_file      = fullfile(mri_dir, 'filename.nii');
        case '0005'
            mri_file      = fullfile(mri_dir, 'filename.nii');
        case '0006'
            mri_file      = fullfile(mri_dir, 'filename.nii');
            % name for bias corrected file:
            mri_file_resl = fullfile(mri_dir, 'specialfile.nii');
        case '0007'
            mri_file      = fullfile(mri_dir, 'filename.nii');
        case '0008'
            mri_file      = fullfile(mri_dir, 'filename.nii');
        case '0009'
            mri_file      = fullfile(mri_dir, 'filename.nii');
        case '0010'
            mri_file      = fullfile(mri_dir, 'filename.nii');
    end

    if strcmp(subj_id, 's0006')
        update_resl = 1;
    else
        update_resl = 0;
        resl_new = 0;
    end

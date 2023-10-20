function EEG_out = detect_bad_samples(EEG,sub)
    
    % Detect samples marked by Rest1 pipeline as bad segments
    EEG_out = EEG;
    bad_samples = find(EEG_out.preproc.rej_segments_mask == 1);
    bad_samples(2,:) = bad_samples(1,:)/250;
    
    % Associate the bad samples to the corresponding task trials
    for t = 1:length(EEG_out.task.PI)
        range = [EEG_out.task.PI(t,1),EEG_out.task.PI(t,1)+EEG_out.task.PI(t,2)];
        bad_samples_bol = bad_samples(2,:) > range(1) & bad_samples(2,:) < range(2);
        bad_samples_idx = bad_samples(1,bad_samples_bol);
        EEG_out.task.PI_badsamples{t} = bad_samples_idx;
    end
    
    for t = 1:length(EEG_out.task.PR)
        range = [EEG_out.task.PR(t,1),EEG_out.task.PR(t,1)+EEG_out.task.PR(t,2)];
        bad_samples_bol = bad_samples(2,:) > range(1) & bad_samples(2,:) < range(2);
        bad_samples_idx = bad_samples(1,bad_samples_bol);
        EEG_out.task.PR_badsamples{t} = bad_samples_idx;
    end
    
    % Add the bad trials detected visually
    EEG_out = manual_outliers(EEG_out,sub);

end


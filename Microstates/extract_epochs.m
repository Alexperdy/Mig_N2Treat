function [EEG_PI_epochs, EEG_PR_epochs] = extract_epochs(EEG)

    % Extract epochs of Pain Imagining (PI)
    EEG_PI_epochs = pop_epoch( EEG, {  'Pain'  }, [0  20], 'newname', 'EEG epochs', 'epochinfo', 'yes');
    EEG_PI_epochs = eeg_checkset( EEG_PI_epochs );
    
    % Extract epochs of Pain Relief (PR) + Last epoch has < than 20 s.
    EEG_PR_epochs = pop_epoch( EEG, {  'Relief'  }, [0  20], 'newname', 'EEG epochs', 'epochinfo', 'yes');
    EEG_PR_lastepoch = pop_epoch( EEG, {  'Relief'  }, [0  EEG.xmax-EEG.task.PR(end,1)], 'newname', 'EEG epochs', 'epochinfo', 'yes');
    EEG_PR_epochs.data(:,:,8) = [EEG_PR_lastepoch.data(:,:,8) zeros(EEG.nbchan,length(EEG_PR_epochs.data)-length(EEG_PR_lastepoch.data))];
    EEG_PR_epochs.trials = EEG_PR_epochs.trials+1;
    EEG_PR_epochs.epoch(8) = EEG_PR_lastepoch.epoch(8);
    EEG_PR_epochs = eeg_checkset( EEG_PR_epochs);

end

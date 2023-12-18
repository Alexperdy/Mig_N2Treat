function EEGout = add_event(EEG, type, latency, duration)

    % Update urevents to match events (to avoid error when adding new events)
    EEGout = EEG;
    EEGout.urevent = EEGout.event;
    EEGout.urevent = rmfield(EEGout.urevent,"urevent");
    
    for j = 1:length(latency)
        
        % Add new event at desired latency (s)
        idx = find(extractfield(EEGout.event,'latency')<=EEGout.srate*latency(j)+1);
        idx = idx(end);
        EEGout = pop_editeventvals(EEGout,'append',{idx [] [] [] []},'changefield',{idx+1 'type' type},'changefield',{idx+1 'latency' latency(j)},'changefield',{idx+1 'duration' duration(j)});
        
        % Update urevents again
        temp = EEGout.event(idx+1).urevent;
        EEGout.event(idx+1).urevent = idx+1;
    
        for i = idx+2:temp
            EEGout.event(i).urevent =  EEGout.event(i).urevent+1;
        end

    end

    EEGout.urevent = EEGout.event;
    EEGout.urevent = rmfield(EEGout.urevent,"urevent");
    EEGout = eeg_checkset( EEGout );

end

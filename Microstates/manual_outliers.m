function EEG_out = manual_outliers(EEG,sub)

EEG_out = EEG;

switch sub
    % controls
    case 19
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 20
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 1 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 25
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 1 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 26
        EEG_out.task.PI_badsamples(2,:) = {1 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 27
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 28
        EEG_out.task.PI_badsamples(2,:) = {0 0 1 0 1 0 0 1};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 29
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 1 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 30
        EEG_out.task.PI_badsamples(2,:) = {1 0 1 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 1 0};

    case 31
        EEG_out.task.PI_badsamples(2,:) = {1 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 33
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 44
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 1 0};

    case 46
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0};

    case 48
        EEG_out.task.PI_badsamples(2,:) = {0 1 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 49
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 51
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    % patients
    case 2
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 3
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 1 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0};

    case 5
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 1 0 1 1 1};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 1 1 1 1 1};

    case 6
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 7
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 1 0 1 1 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 1 0 0 0 0 0};

    case 8
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 9
        EEG_out.task.PI_badsamples(2,:) = {1 0 1 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 1 1 0 1};

    case 12
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 1 0};

    case 13
        EEG_out.task.PI_badsamples(2,:) = {0 1 0 1 1 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 1 0 1};

    case 34
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 38
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 0 0};

    case 41
        EEG_out.task.PI_badsamples(2,:) = {0 0 0 0 0 1 1 0};
        EEG_out.task.PR_badsamples(2,:) = {0 0 0 0 0 0 1 0};

    case 43
        EEG_out.task.PI_badsamples(2,:) = {0 1 0 0 0 1 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 1 0 0 0 1 0 0};

    case 45
        EEG_out.task.PI_badsamples(2,:) = {0 1 1 0 1 0 0 0};
        EEG_out.task.PR_badsamples(2,:) = {0 1 1 1 1 0 1 0};
end


end
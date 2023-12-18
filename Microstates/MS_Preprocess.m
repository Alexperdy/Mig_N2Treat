% % input: preprocessed EEG from Rest1 (example: sub-control019_ses-midcycle_task-mentalizing_eeg_pre-Rest1)

%% Set paths
clearvars; clc; close all

% if patient: load paths and variables for file names of patients
[setpath,setfil,task,type,ses]=load_paths('patients');
sub = [2,3,5,6,7,8,9,12,13,34,38,41,43,45];

% If control: load paths and variables for file names of controls
[setpath,setfil,task,type,ses]=load_paths('controls');
sub = [19,20,25,26,27,28,29,30,31,33,44,46,48,49,51];


%% Remove/inteprolate outliers and save the data
% Uses capmean: interpolate bad samples (mean +/- 4 std) with window = 2;
 for i = sub

    % Load EEG data
    EEG_out = pop_loadset(sprintf(setfil,type,i,ses,task),sprintf(setpath,type,i));
    EEG_outlier = run_outlier_rej(EEG_out,'capmean',2);
    pop_saveset(EEG_outlier, strcat(sprintf(setfil,type,i,ses,task),'_outliers'), sprintf(setpath,type,i));

    clear EEG_out EEG_outlier
 end

%% Load data for sub(i), add task timings & extract epochs:

for i = sub
    % Load EEG data
    EEG_sub = pop_loadset(strcat(sprintf(setfil,type,i,ses,task),'_outliers.set'),sprintf(setpath,type,i));

    % Add task timings: (Pain Imag 20 s -> Pain Relief 20 s)x8

    task_file = strcat(sprintf(setpath,type,i),sprintf('\\sub-%s%03d_ses-%s_task-mentalizing_acq-ep2d_p2_s3_stim-pain_correct.txt',type,i,ses));
    task_matrix = readmatrix(task_file);
    PI_lat = task_matrix(:,1)';
    PI_dur = task_matrix(:,2)';
    PR_lat= PI_lat+PI_dur;
    PR_dur = [PI_lat(2:end)-PR_lat(1:end-1),EEG_sub.xmax-PR_lat(end)];
    EEG_sub.task.PI = [PI_lat',PI_dur'];
    EEG_sub.task.PR = [PR_lat',PR_dur'];

    EEG_task = add_event(EEG_sub,'Pain',PI_lat,PI_dur);
    EEG_task = add_event(EEG_task,'Relief',PR_lat,PR_dur);
    
    % Detect bad samples:
    EEG_task = detect_bad_samples(EEG_task,i);

    % Extract epochs of Pain Imagining (PI) & Pain Relief (PR)
    [EEG_PI_epochs, EEG_PR_epochs] = extract_epochs(EEG_task);

    % Save EEG data with task events & extracted epochs
    pop_saveset(EEG_task, strcat(sprintf(setfil,type,i,ses,task),'_task'), sprintf(setpath,type,i));
    pop_saveset(EEG_PI_epochs, strcat(sprintf(setfil,type,i,ses,task),'_epochPI'), sprintf(setpath,type,i));
    pop_saveset(EEG_PR_epochs, strcat(sprintf(setfil,type,i,ses,task),'_epochPR'), sprintf(setpath,type,i));
    clear EEG_sub EEG_task EEG_PI_epochs EEG_PR_epochs PI_lat PI_dur PR_lat PR_dur task_matrix task_file
end
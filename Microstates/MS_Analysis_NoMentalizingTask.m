%% Define the basic parameters

clear all
close all
clc

LowCutFilter  =  2;
HighCutFilter = 20;
FilterCoefs   = 2000;


% EEG inside sanner (resting state)
%SavePath = 'D:\\Universidade\\5º Ano 1º Semestre\\Thesis\\MATLAB\\DATA\\MS_DATA\\Resting_state\\EEG_data\\sub-%s%03d\\ses-%s\\preprocessed\\sub-%s%03d_ses-%s_task-rest_eeg_preproc';
%SavePath2 = 'D:\\Universidade\\5º Ano 1º Semestre\\Thesis\\MATLAB\\DATA\\MS_DATA\\Resting_state\\';

% EEG outside sanner: calibration task
SavePath = 'D:\\Universidade\\5º Ano 1º Semestre\\Thesis\\MATLAB\\DATA\\MS_DATA\\Calibration_task\\Outside_scanner\\eeg_task-eegcalibrationout\\sub-%s%03d\\ses-%s\\preprocessed\\sub-%s%03d_ses-%s_task-eegcalibrationout_eeg_preproc';
SavePath2 = 'D:\\Universidade\\5º Ano 1º Semestre\\Thesis\\MATLAB\\DATA\\MS_DATA\\Calibration_task\\\Outside_scanner\\';
EventPath = 'D:\\Universidade\\5º Ano 1º Semestre\\Thesis\\MATLAB\\DATA\\MS_DATA\\Calibration_task\\Outside_scanner\\eeg_task-eegcalibrationout\\sub-%s%03d\\ses-%s\\preprocessed\\events_sub-%s%03d_ses-%s_task-eegcalibrationout';

% EEG inside sanner: calibration task
%SavePath = 'D:\\Universidade\\5º Ano 1º Semestre\\Thesis\\MATLAB\\DATA\\MS_DATA\\Calibration_task\\Inside_scanner\\eeg_task-eegcalibration\\sub-%s%03d\\ses-%s\\preprocessed\\sub-%s%03d_ses-%s_task-eegcalibration_eeg_preproc';
%SavePath2 = 'D:\\Universidade\\5º Ano 1º Semestre\\Thesis\\MATLAB\\DATA\\MS_DATA\\Calibration_task\\Inside_scanner\\';


nGroups = 2;
%% Read the data

eeglab

AllSubjects = [];

%[setpath,setfil,task,type,ses]=load_paths('controls');
controls = [19,20,25,26,27,28,29,30,31,33,44,46,48,49,51];
GroupIndex{1} = [];

for i = controls
    tmpEEG = load(sprintf(SavePath,'control',i,'midcycle','control',i,'midcycle'));
    tmpEEG = tmpEEG.EEG;
    tmpEEG = pop_select(tmpEEG,'rmchannel',32);
    % KEEP THIS SECTION ONLY WHEN USING EYES OPEN DATA
    events = load(sprintf(EventPath,'control',i,'midcycle','control',i,'midcycle'));
    events = events.events;
    eventsEO = events.EO(1:31,:);
    tmpEEG.data = eventsEO;
    %

    [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'gui','off'); % And make this a new set
    %tmpEEG = pop_clean_rawdata(tmpEEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',5,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
    tmpEEG = pop_reref(tmpEEG, []); % Make things average reference
    tmpEEG = pop_eegfiltnew(tmpEEG, LowCutFilter,HighCutFilter, FilterCoefs, 0, [], 0); % bandpass-filter 2-20Hz
    tmpEEG.group = 'Group_Control'; % Controls = Group 1
    tmpEEG.setname = sprintf('sub-control%03d_ses-midcycle_task-rest_eeg',i);
    [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, tmpEEG, CURRENTSET); % Store the thing
    GroupIndex{1} = [GroupIndex{1} CURRENTSET]; % And keep track of the group
    AllSubjects = [AllSubjects CURRENTSET];

end

%[setpath,setfil,task,type,ses]=load_paths('patients');
patients = [2,3,5,6,7,8,9,12,13,34,38,41,43,45];
GroupIndex{2} = [];

for i = patients
    tmpEEG = load(sprintf(SavePath,'patient',i,'interictal','patient',i,'interictal'));
    tmpEEG = tmpEEG.EEG;
    tmpEEG = pop_select(tmpEEG,'rmchannel',32);
    % KEEP THIS SECTION ONLY WHEN USING EYES OPEN DATA
    events = load(sprintf(EventPath,'patient',i,'interictal','patient',i,'interictal'));
    events = events.events;
    eventsEO = events.EO(1:31,:);
    tmpEEG.data = eventsEO;
    %
    [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'gui','off'); % And make this a new set
    %tmpEEG = pop_clean_rawdata(tmpEEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',5,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
    tmpEEG = pop_reref(tmpEEG, []); % Make things average reference
    tmpEEG = pop_eegfiltnew(tmpEEG, LowCutFilter,HighCutFilter, FilterCoefs, 0, [], 0); % bandpass-filter 2-20Hz
    tmpEEG.group = 'Group_Patients'; % Patients = Group 2
    tmpEEG.setname = sprintf('sub-patient%03d_ses-midcycle_task-rest_eeg',i);
    [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, tmpEEG, CURRENTSET); % Store the thing
    GroupIndex{2} = [GroupIndex{2} CURRENTSET]; % And keep track of the group
    AllSubjects = [AllSubjects CURRENTSET];
end

eeglab redraw

%% Identify microstates for each individual

% Define the parameters for clustering: 4 Classes, Uses only GFP Peaks,
% uses all GFP peaks, and uses AAHC algorithm
ClustPars = struct('MinClasses',4,'MaxClasses',4,'GFPPeaks',1,'IgnorePolarity',1,'MaxMaps',Inf,'Restarts',50', 'UseAAHC',1,'Normalize',0);

% Loop across all subjects to identify the individual clusters
for i = 1:numel(AllSubjects) 
    EEG = eeg_retrieve(ALLEEG,AllSubjects(i)); % the EEG we want to work with
    fprintf(1,'Clustering dataset %s (%i/%i)\n',EEG.setname,i,numel(AllSubjects )); % Some info for the impatient user
    EEG = pop_FindMSTemplates(EEG, ClustPars); % This is the actual clustering within subjects
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, AllSubjects(i)); % Done, we just need to store this
end

eeglab redraw   

%% Now we combine the microstate maps across subjects and sort the mean

% First, we load a set of normative maps to orient us later

templatepath = fullfile(fileparts(which('eegplugin_Microstates.m')),'Templates');

EEG = pop_loadset('filename','Normative microstate template maps Neuroimage 2002.set','filepath',templatepath);
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); % And make this a new set

% And we have a look at it
NormativeTemplateIndex = CURRENTSET;
pop_ShowIndMSMaps(ALLEEG(NormativeTemplateIndex), 4); 
drawnow;

% Now we go into averaging within each group
for Group = 1:nGroups
    % The mean of group X
    EEG = pop_CombMSTemplates(ALLEEG, GroupIndex{Group}, 0, 0, sprintf('GrandMean Group %i',Group));
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, numel(ALLEEG)+1,'gui','off'); % Make a new set
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % and store it
    GrandMeanIndex(Group) = CURRENTSET; % And keep track of it
end

% % We automatically sort the means based on a template from the literature
for Group = 1:nGroups
    [ALLEEG,EEG] = pop_SortMSTemplates(ALLEEG, GrandMeanIndex(Group), 1, NormativeTemplateIndex);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, GrandMeanIndex(Group));
end

eeglab redraw

%% Save data
eeglab redraw

save(strcat(SavePath2,'ALLEEG.mat'),'ALLEEG')
save(strcat(SavePath2,'AllSubjects.mat'),'AllSubjects')
save(strcat(SavePath2,'GrandMeanIndex.mat'),'GrandMeanIndex')
save(strcat(SavePath2,'GroupIndex.mat'),'GroupIndex')
save(strcat(SavePath2,'NormativeTemplateIndex.mat'),'NormativeTemplateIndex')

%% Load data

eeglab
load(strcat(SavePath2,'ALLEEG.mat'))
load(strcat(SavePath2,'AllSubjects.mat'))
load(strcat(SavePath2,'GrandMeanIndex.mat'))
load(strcat(SavePath2,'GroupIndex.mat'))
load(strcat(SavePath2,'NormativeTemplateIndex.mat'))

eeglab redraw

%% Compute grand-grand mean

% Now we want the grand-grand mean, based on the group means, if there is
% more than one group
if nGroups > 1
    EEG = pop_CombMSTemplates(ALLEEG, GrandMeanIndex, 1, 0, 'GrandGrandMean');
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, numel(ALLEEG)+1,'gui','off'); % Make a new set
    GrandGrandMeanIndex = CURRENTSET; % and keep track of it
else
    GrandGrandMeanIndex = GrandMeanIndex(1);
end

% We automatically sort the grandgrandmean based on a template from the literature
[ALLEEG,EEG] = pop_SortMSTemplates(ALLEEG, GrandGrandMeanIndex, 1, NormativeTemplateIndex);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, GrandGrandMeanIndex);

eeglab redraw
%% And we sort things out over means and subjects
% Now, that we have mean maps, we use them to sort the individual templates
nGroups = 2;

% First, the sequence of the two group means has be adjusted based on the
% grand grand mean

if nGroups > 1
    ALLEEG = pop_SortMSTemplates(ALLEEG, GrandMeanIndex, 1, GrandGrandMeanIndex);
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % and store it
end

% Then, we sort the individuals based on their group means
for Group = 1:nGroups
    ALLEEG = pop_SortMSTemplates(ALLEEG, GroupIndex{Group}, 0, GrandMeanIndex(Group)); % Group 1
    [ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % and store it
end

eeglab redraw

%% Stats part
%SavePath = 'D:\\Universidade\\5º Ano 1º Semestre\\Thesis\\MATLAB\\MS_DATA\\AAHC_2_20_Inf\\';
nGroups = 2;
% These are the paramters for the fitting based on GFP peaks only
FitPars = struct('nClasses',4,'lambda',1,'b',25,'PeakFit',1, 'BControl',1,'Rectify',0,'Normalize',0);

% Using the individual templates
pop_QuantMSTemplates(ALLEEG, AllSubjects, 0, FitPars, []                   , fullfile(SavePath2,'ResultsFromIndividualTemplates_PeakFit2.xlsx'));

%% Eventually export the individual microstate maps to do statistics in Ragu
pop_RaguMSTemplates(ALLEEG, AllSubjects);
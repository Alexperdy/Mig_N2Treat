% More info: http://www.thomaskoenig.ch/index.php/work/software/microstates-in-eeglab/getting-started
% More info for functions: help function_name

% Before running this code:
% run the script MS_preprocess to obtain the preprocessed eeg data with outliers 
% interpolated (files eegdata_outliers) and eeg data with pain imagery (PI) 
% and pain relief (PR) epochs separated (files eegdata_epochPI and eegdata_epochPR)

%% Define the basic parameters

clear all
close all
clc

% Band-pass filter (2-20Hz)
LowCutFilter  =  2;
HighCutFilter = 20;
FilterCoefs   = 2000;

% Path to save data
SavePath = 'D:\\Universidade\\5º Ano 1º Semestre\\Thesis\\MATLAB\\DATA\\MS_DATA\\AAHC\\';

% Path with eeg data: eeg already preprocessed with outliers interpolated
eeg_path = 'D:\\Universidade\\5º Ano 1º Semestre\\Thesis\\MATLAB\\EEG_data\\sub-%s%03d\\sub-%s%03d_ses-%s_task-mentalizing_eeg_pre-Rest1_outliers.set';

% Number of groups: 2 (Control & Patients)
nGroups = 2;


%% Read the data

eeglab

AllSubjects = [];

% load paths and variables for file names of controls
[setpath,setfil,task,type,ses]=load_paths('controls');
controls = [19,20,25,26,27,28,29,30,31,33,44,46,48,49,51];
GroupIndex{1} = [];

for i = controls
    tmpEEG = pop_loadset(sprintf(eeg_path,'control',i,'control',i,'midcycle'));
    [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'gui','off'); % And make this a new set
    tmpEEG = pop_clean_rawdata(tmpEEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',5,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
    tmpEEG = pop_reref(tmpEEG, []); % Make things average reference
    tmpEEG = pop_eegfiltnew(tmpEEG, LowCutFilter,HighCutFilter, FilterCoefs, 0, [], 0); % bandpass-filter 2-20Hz
    tmpEEG.group = 'Group_Control'; % Controls = Group 1
    tmpEEG.setname = sprintf('sub-control%03d_ses-midcycle_task-mentalizing_eeg',i);
    [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, tmpEEG, CURRENTSET); % Store the data
    GroupIndex{1} = [GroupIndex{1} CURRENTSET]; % And keep track of the group
    AllSubjects = [AllSubjects CURRENTSET];

end

% load paths and variables for file names of patients
[setpath,setfil,task,type,ses]=load_paths('patients');
patients = [2,3,5,6,7,8,9,12,13,34,38,41,43,45];
GroupIndex{2} = [];

for i = patients
    tmpEEG = pop_loadset(sprintf(eeg_path,'patient',i,'patient',i,'interictal'));
    [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'gui','off'); % And make this a new set
    tmpEEG = pop_clean_rawdata(tmpEEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',5,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
    tmpEEG = pop_reref(tmpEEG, []); % Make things average reference
    tmpEEG = pop_eegfiltnew(tmpEEG, LowCutFilter,HighCutFilter, FilterCoefs, 0, [], 0); % bandpass-filter 2-20Hz
    tmpEEG.group = 'Group_Patients'; % Patients = Group 2
    tmpEEG.setname = sprintf('sub-patient%03d_ses-midcycle_task-mentalizing_eeg',i);
    [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, tmpEEG, CURRENTSET); % Store the data
    GroupIndex{2} = [GroupIndex{2} CURRENTSET]; % And keep track of the group
    AllSubjects = [AllSubjects CURRENTSET];
end

eeglab redraw

%% Identify microstates for each individual

% Define the parameters for clustering
ClustPars = struct('MinClasses',4,'MaxClasses',4,'GFPPeaks',1,'IgnorePolarity',1,'MaxMaps',Inf,'Restarts',50', 'UseAAHC',1,'Normalize',0);

% Loop across all subjects to identify the individual clusters
for i = 1:numel(AllSubjects ) 
    EEG = eeg_retrieve(ALLEEG,AllSubjects(i)); % the EEG we want to work with
    fprintf(1,'Clustering dataset %s (%i/%i)\n',EEG.setname,i,numel(AllSubjects )); 
    EEG = pop_FindMSTemplates(EEG, ClustPars); % Clustering within each subject
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, AllSubjects(i)); % Store the data
end

eeglab redraw   

%% Combine the microstate maps across subjects and sort the mean

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

%% Interactive sorting of the mean maps - Run one by one

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',GrandMeanIndex(1),'study',0);
[ALLEEG,EEG] = pop_ShowIndMSMaps(EEG, 4, 1, ALLEEG);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, GrandMeanIndex(1));

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',GrandMeanIndex(2),'study',0);
[ALLEEG,EEG] = pop_ShowIndMSMaps(EEG, 4, 1, ALLEEG);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, GrandMeanIndex(2));

eeglab redraw

%% Save data
eeglab redraw

save(strcat(SavePath,'ALLEEG.mat'),'ALLEEG')
save(strcat(SavePath,'AllSubjects.mat'),'AllSubjects')
save(strcat(SavePath,'GrandMeanIndex.mat'),'GrandMeanIndex')
save(strcat(SavePath,'GroupIndex.mat'),'GroupIndex')
save(strcat(SavePath,'NormativeTemplateIndex.mat'),'NormativeTemplateIndex')

%% Load data

eeglab
load(strcat(SavePath,'ALLEEG.mat'))
load(strcat(SavePath,'AllSubjects.mat'))
load(strcat(SavePath,'GrandMeanIndex.mat'))
load(strcat(SavePath,'GroupIndex.mat'))
load(strcat(SavePath,'NormativeTemplateIndex.mat'))

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load and preprocess the data separated by PR and PI periods 
% We can want to fit the obtained maps obtained by using all
% the data (what we did above!) to the PR and PI data.

LowCutFilter  =  2;
HighCutFilter = 20;
FilterCoefs   = 2000;

idxStats = [];

[setpath,setfil,task,type,ses]=load_paths('controls');
controls = [19,20,25,26,27,28,29,30,31,33,44,46,48,49,51];

for c = 1:numel(controls)
    % Add PI
    i = controls(c);
    tmpEEG = pop_loadset(strcat(sprintf(setfil,type,i,ses,task),'_epochPI.set'),sprintf(setpath,type,i));
    [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'gui','off'); % And make this a new set
    tmpEEG = pop_select(tmpEEG,'trial',find(horzcat(tmpEEG.task.PI_badsamples{2,:})==0)); % remove bad trials
    tmpEEG = pop_reref(tmpEEG, []); % Make things average reference
    tmpEEG = pop_eegfiltnew(tmpEEG, LowCutFilter,HighCutFilter, FilterCoefs, 0, [], 0); % bandpass-filter 2-20Hz
    tmpEEG.msinfo = ALLEEG(c).msinfo;
    tmpEEG.group = 'Group_Control_PI'; % Controls = Group 1
    tmpEEG.setname = sprintf('sub-control%03d_ses-midcycle_task-mentalizing_eeg_PI',i);
    [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, tmpEEG, CURRENTSET); % Store the thing
    idxStats = [idxStats CURRENTSET];

    % Add PR
    tmpEEG = pop_loadset(strcat(sprintf(setfil,type,i,ses,task),'_epochPR.set'),sprintf(setpath,type,i));
    [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'gui','off'); % And make this a new set
    tmpEEG = pop_select(tmpEEG,'trial',find(horzcat(tmpEEG.task.PR_badsamples{2,:})==0));
    tmpEEG = pop_reref(tmpEEG, []); % Make things average reference
    tmpEEG = pop_eegfiltnew(tmpEEG, LowCutFilter,HighCutFilter, FilterCoefs, 0, [], 0); % bandpass-filter 2-20Hz
    tmpEEG.msinfo = ALLEEG(c).msinfo;
    tmpEEG.group = 'Group_Control_PR'; % Controls = Group 1
    tmpEEG.setname = sprintf('sub-control%03d_ses-midcycle_task-mentalizing_eeg_PR',i);
    [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, tmpEEG, CURRENTSET); % Store the thing
    idxStats = [idxStats CURRENTSET];
end

[setpath,setfil,task,type,ses]=load_paths('patients');
patients = [2,3,5,6,7,8,9,12,13,34,38,41,43,45];

for p = 1:numel(patients)
    % Add PI
    i = patients(p);
    tmpEEG = pop_loadset(strcat(sprintf(setfil,type,i,ses,task),'_epochPI.set'),sprintf(setpath,type,i));
    [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'gui','off'); % And make this a new set
    tmpEEG = pop_select(tmpEEG,'trial',find(horzcat(tmpEEG.task.PI_badsamples{2,:})==0));
    tmpEEG = pop_reref(tmpEEG, []); % Make things average reference
    tmpEEG = pop_eegfiltnew(tmpEEG, LowCutFilter,HighCutFilter, FilterCoefs, 0, [], 0); % bandpass-filter 2-20Hz
    tmpEEG.msinfo = ALLEEG(15+p).msinfo;
    tmpEEG.group = 'Group_Patients_PI'; % Patients = Group 2
    tmpEEG.setname = sprintf('sub-patient%03d_ses-midcycle_task-mentalizing_eeg_PI',i);
    [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, tmpEEG, CURRENTSET); % Store the thing
    idxStats = [idxStats CURRENTSET];

    % Add PR
    tmpEEG = pop_loadset(strcat(sprintf(setfil,type,i,ses,task),'_epochPR.set'),sprintf(setpath,type,i));
    [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'gui','off'); % And make this a new set
    tmpEEG = pop_select(tmpEEG,'trial',find(horzcat(tmpEEG.task.PR_badsamples{2,:})==0));
    tmpEEG = pop_reref(tmpEEG, []); % Make things average reference
    tmpEEG = pop_eegfiltnew(tmpEEG, LowCutFilter,HighCutFilter, FilterCoefs, 0, [], 0); % bandpass-filter 2-20Hz
    tmpEEG.msinfo = ALLEEG(15+p).msinfo;
    tmpEEG.group = 'Group_Patients_PR'; % Controls = Group 1
    tmpEEG.setname = sprintf('sub-patient%03d_ses-midcycle_task-mentalizing_eeg_PR',i);
    [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, tmpEEG, CURRENTSET); % Store the thing
    idxStats = [idxStats CURRENTSET];
end

eeglab redraw

%% Stats part - Fit the data and obtain microstate parameters
%SavePath = 'D:\\Universidade\\5º Ano 1º Semestre\\Thesis\\MATLAB\\MS_DATA\\AAHC_2_20_Inf\\';
nGroups = 2;
% These are the paramters for the fitting based on GFP peaks only
FitPars = struct('nClasses',4,'lambda',1,'b',25,'PeakFit',1, 'BControl',1,'Rectify',0,'Normalize',0);

% Using the individual templates
pop_QuantMSTemplates(ALLEEG, idxStats, 0, FitPars, []                   , fullfile(SavePath,'ResultsFromIndividualTemplates_PeakFit_30_10.xlsx'));

%% Below is the code to fit the data using the group MS, or the grand mean maps, 
%% or the literature maps

%% Using group mean templates

%GroupIndexStats{1} = idxStats(1:15*2); 
%GroupIndexStats{2} = idxStats(15*2+1:end);

%for Group = 1:nGroups
%    pop_QuantMSTemplates(ALLEEG, GroupIndexStats{Group}, 1, FitPars, GrandMeanIndex(Group), fullfile(SavePath,sprintf('ResultsFromMeanTemplateGroup_%i_AllData_30_10.xlsx',Group)));
%end

%% And using the grand grand mean template
% pop_QuantMSTemplates(ALLEEG, idxStats, 1, FitPars, GrandGrandMeanIndex, fullfile(SavePath,'ResultsFromGrandGrandMeanTemplate_AllData_30_10.xlsx'));

%% And finally, based on the normative maps from 2002
% pop_QuantMSTemplates(ALLEEG, idxStats, 1, FitPars, NormativeTemplateIndex, fullfile(SavePath,'ResultsFromNormativeTemplate2002.xlsx'));

%% Eventually export the individual microstate maps to do statistics in Ragu
pop_RaguMSTemplates(ALLEEG, AllSubjects);
%% Define the basic parameters

clear all
close all
clc

LowCutFilter  =  2;
HighCutFilter = 20;
FilterCoefs   = 2000;

% Path where the data is saved
SavePath = 'D:\\Universidade\\5º Ano 1º Semestre\\Thesis\\MATLAB\\DATA\\MS_DATA\\AAHC\\';

% Number of groups: 2 (Control & Patients)
nGroups = 2;

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

%% Add the original data and assign them the corresponding MS Maps
LowCutFilter  =  2;
HighCutFilter = 20;
FilterCoefs   = 2000;

idxStats = [];

[setpath,setfil,task,type,ses]=load_paths('controls');
controls = [19,20,25,26,27,28,29,30,31,33,44,46,48,49,51];

for c = 1:numel(controls)

    i = controls(c);
    tmpEEG = pop_loadset(strcat(sprintf(setfil,type,i,ses,task),'_outliers.set'),sprintf(setpath,type,i));
    [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'gui','off'); % And make this a new set
    tmpEEG = pop_reref(tmpEEG, []); % Make things average reference
    tmpEEG = pop_eegfiltnew(tmpEEG, LowCutFilter,HighCutFilter, FilterCoefs, 0, [], 0); % bandpass-filter 2-20Hz
    tmpEEG.msinfo = ALLEEG(c).msinfo;
    tmpEEG.group = 'Group_Control'; % Controls = Group 1
    tmpEEG.setname = sprintf('sub-control%03d',i);
    [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, tmpEEG, CURRENTSET); % Store the thing
    idxStats = [idxStats CURRENTSET];
end

[setpath,setfil,task,type,ses]=load_paths('patients');
patients = [2,3,5,6,7,8,9,12,13,34,38,41,43,45];

for p = 1:numel(patients)
    % Add PI
    i = patients(p);
    tmpEEG = pop_loadset(strcat(sprintf(setfil,type,i,ses,task),'_outliers.set'),sprintf(setpath,type,i));
    [ALLEEG, tmpEEG, CURRENTSET] = pop_newset(ALLEEG, tmpEEG, 0,'gui','off'); % And make this a new set
    tmpEEG = pop_reref(tmpEEG, []); % Make things average reference
    tmpEEG = pop_eegfiltnew(tmpEEG, LowCutFilter,HighCutFilter, FilterCoefs, 0, [], 0); % bandpass-filter 2-20Hz
    tmpEEG.msinfo = ALLEEG(15+p).msinfo;
    tmpEEG.group = 'Group_Patients'; % Patients = Group 2
    tmpEEG.setname = sprintf('sub-patient%03d',i);
    [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, tmpEEG, CURRENTSET); % Store the thing
    idxStats = [idxStats CURRENTSET];
end

eeglab redraw

%% Fit MS Maps and obtain the time series of each MS

% Fit: Fit only the GFP Peaks and interpolate
FitPars = struct('nClasses',4,'lambda',1,'b',25,'PeakFit',1, 'BControl',1,'Rectify',0,'Normalize',0);

for i = 1:numel(idxStats)
    TheEEG = eeg_retrieve(ALLEEG,idxStats(i));
    Maps = TheEEG.msinfo.MSMaps(4).Maps;
    [MSClass,gfp,fit] = AssignMStates(TheEEG,Maps,FitPars,TheEEG.msinfo.ClustPar.IgnorePolarity);
    TheEEG.msinfo.FitInfo.MSClass = MSClass;
    TheEEG.msinfo.FitInfo.gfp = gfp;
    TheEEG.msinfo.FitInfo.expvar = fit;

    TheEEG.msinfo.FitInfo.MSTimeSeries = zeros(4,length(MSClass));
    TheEEG.msinfo.FitInfo.MSTimeSeries(1,find(MSClass == 1)) = 1;
    TheEEG.msinfo.FitInfo.MSTimeSeries(2,find(MSClass == 2)) = 1;
    TheEEG.msinfo.FitInfo.MSTimeSeries(3,find(MSClass == 3)) = 1;
    TheEEG.msinfo.FitInfo.MSTimeSeries(4,find(MSClass == 4)) = 1;

    TheEEG.msinfo.FitInfo.MSCorrSeries = zeros(4,length(MSClass));
    TheEEG.msinfo.FitInfo.MSCorrSeries(1,:) = abs(corr(Maps(1,:)',TheEEG.data));
    TheEEG.msinfo.FitInfo.MSCorrSeries(2,:) = abs(corr(Maps(2,:)',TheEEG.data));
    TheEEG.msinfo.FitInfo.MSCorrSeries(3,:) = abs(corr(Maps(3,:)',TheEEG.data));
    TheEEG.msinfo.FitInfo.MSCorrSeries(4,:) = abs(corr(Maps(4,:)',TheEEG.data));

    mkdir(strcat(SavePath,'PeakFit\\EEG data\\',TheEEG.setname));
    save(strcat(SavePath,'PeakFit\\EEG data\\',TheEEG.setname,'\\',TheEEG.setname,'_MS'),'TheEEG');
end




%% PLOTS

t = TheEEG.times/1000;
figure()
subplot(4,1,1)
plot(t,TheEEG.msinfo.FitInfo.MSTimeSeries(1,:),'LineWidth',2)
ylim([0,1.5])
xlim([100 105])
legend('Microstate A','FontSize',14,'FontWeight','bold')
subplot(4,1,2)
plot(t,TheEEG.msinfo.FitInfo.MSTimeSeries(2,:),'r','LineWidth',2)
ylim([0,1.5])
xlim([100 105])
legend('Microstate B','FontSize',14,'FontWeight','bold')
subplot(4,1,3)
plot(t,TheEEG.msinfo.FitInfo.MSTimeSeries(3,:),'g','LineWidth',2)
ylim([0,1.5])
xlim([100 105])
legend('Microstate C','FontSize',14,'FontWeight','bold')
subplot(4,1,4)
plot(t,TheEEG.msinfo.FitInfo.MSTimeSeries(4,:),'m','LineWidth',2)
ylim([0,1.5])
xlim([100 105])
legend('Microstate D','FontSize',14,'FontWeight','bold')
xlabel('Time (s)','FontSize',14,'FontWeight','bold')
sgtitle('Occurrence Time Series')

t = TheEEG.times/1000;
figure()
subplot(4,1,1)
plot(t,TheEEG.msinfo.FitInfo.MSCorrSeries(1,:),'LineWidth',2)
% hold on
% plot(t,TheEEG.msinfo.FitInfo.MSTimeSeries(1,:),'k','LineWidth',2)
ylim([0,1])
xlim([100 105])
ylabel('Spatial Correlation')
legend('Microstate A','FontSize',14,'FontWeight','bold')
subplot(4,1,2)
plot(t,TheEEG.msinfo.FitInfo.MSCorrSeries(2,:),'r','LineWidth',2)
% hold on
% plot(t,TheEEG.msinfo.FitInfo.MSTimeSeries(2,:),'k','LineWidth',2)
ylim([0,1])
xlim([100 105])
ylabel('Spatial Correlation')
legend('Microstate B','FontSize',14,'FontWeight','bold')
subplot(4,1,3)
plot(t,TheEEG.msinfo.FitInfo.MSCorrSeries(3,:),'g','LineWidth',2)
% hold on
% plot(t,TheEEG.msinfo.FitInfo.MSTimeSeries(3,:),'k','LineWidth',2)
ylim([0,1])
xlim([100 105])
ylabel('Spatial Correlation')
legend('Microstate C','FontSize',14,'FontWeight','bold')
subplot(4,1,4)
plot(t,TheEEG.msinfo.FitInfo.MSCorrSeries(4,:),'m','LineWidth',2)
% hold on
% plot(t,TheEEG.msinfo.FitInfo.MSTimeSeries(4,:),'k','LineWidth',2)
ylim([0,1])
xlim([100 105])
xlabel('Spatial Correlation')
legend('Microstate D','FontSize',14,'FontWeight','bold')
xlabel('Time (s)','FontSize',14,'FontWeight','bold')
ylabel('Spatial Correlation')
sgtitle('Spatial Correlation Time Series')


t = TheEEG.times/1000;
figure()
plot(t,TheEEG.msinfo.FitInfo.MSCorrSeries(1,:),'b','LineWidth',2,'LineStyle',"--")
hold on
plot(t,TheEEG.msinfo.FitInfo.MSCorrSeries(2,:),'r','LineWidth',2)
ylim([0,1])
xlim([100 101])
ylabel('Spatial Correlation','FontSize',14,'FontWeight','bold')
legend('Microstate A','Microstate B','FontSize',14,'FontWeight','bold')
xlabel('Time (s)','FontSize',14,'FontWeight','bold')

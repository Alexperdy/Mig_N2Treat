function [setpath,setfil,task,type,ses] = load_paths(subj)

setpath = 'D:\\Universidade\\5º Ano 1º Semestre\\Thesis\\MATLAB\\EEG_data\\sub-%s%03d';
%setpath = 'C:\\Users\\alexp\\OneDrive\\Ambiente de Trabalho\\Thesis\\MATLAB\\EEG_data\\sub-%s%03d';
setfil = 'sub-%s%03d_ses-%s_task-%s_eeg_pre-Rest1';
task = 'mentalizing';

switch subj
    case 'patients'
        % If patient
        type = 'patient';
        ses = 'interictal';
    case 'controls'
        type = 'control';
        ses = 'midcycle';

end
      
end
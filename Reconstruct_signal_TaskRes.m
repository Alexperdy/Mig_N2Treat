% controls
% path_design = 'D:\\Universidade\\5º Ano 1º Semestre\\FUNC_DATA\\design_matrices\\sub-control%03d\\task-%s';
% path_pe = 'D:\\Universidade\\5º Ano 1º Semestre\\FUNC_DATA\\PEs\\sub-control%03d\\task-%s';
% path_data = 'D:\\Universidade\\5º Ano 1º Semestre\\FUNC_DATA\\filtered_func\\sub-control%03d\\task-%s';
% path_res = 'D:\\Universidade\\5º Ano 1º Semestre\\FUNC_DATA\\residuos\\sub-control%03d\\task-%s';
% path_final = 'D:\\Universidade\\5º Ano 1º Semestre\\FUNC_DATA\\reconstruct_data\\sub-control%03d\\task-%s';
% sub = [19,20,25,26,27,28,29,30,31,33,44,46,48,49,51];

% patients
path_design = 'D:\\Universidade\\5º Ano 1º Semestre\\FUNC_DATA\\design_matrices\\sub-patient%03d\\task-%s';
path_pe = 'D:\\Universidade\\5º Ano 1º Semestre\\FUNC_DATA\\PEs\\sub-patient%03d\\task-%s';
path_data = 'D:\\Universidade\\5º Ano 1º Semestre\\FUNC_DATA\\filtered_func\\sub-patient%03d\\task-%s';
path_res = 'D:\\Universidade\\5º Ano 1º Semestre\\FUNC_DATA\\residuos\\sub-patient%03d\\task-%s';
path_final = 'D:\\Universidade\\5º Ano 1º Semestre\\FUNC_DATA\\reconstruct_data\\sub-patient%03d\\task-%s';
sub = [2,3,5,6,7,8,9,12,13,34,38,41,43,45];


task = 'mentalizing';

%%

for s = 1:length(sub)
    l = sub(s);
    
    % Individualize regressors
    design_mat = dlmread(strcat(sprintf(path_design,l,task), '\design.mat'));
    
    ev1 = design_mat(:,1); % task
    ev2 = design_mat(:,2); % task time derivative

    clear design_mat;

    % Load PEs
    pe1 = niftiread(strcat(sprintf(path_pe,l,task),'\pe1','.nii.gz'));
    pe2 = niftiread(strcat(sprintf(path_pe,l,task),'\pe2','.nii.gz'));

    % Load Residuals
    res = niftiread(strcat(sprintf(path_res,l,task),'\res4d','.nii.gz'));

    % Load header
    header = niftiinfo(strcat(sprintf(path_data,l,task),'\filtered_func_data.nii.gz'));

    % generate model with PEs, EVs and residuals from FEAT
    timepoints = size(res,4);
    model = [];
    for t =1:timepoints
        model(:,:,:,t) = pe1.*ev1(t) + pe2.*ev2(t);
    end

    final_data = model + res;
    niftiwrite(final_data,strcat(sprintf(path_final,l,task),'\final_data_wres'),header);
end
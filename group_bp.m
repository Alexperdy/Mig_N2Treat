function [subjs_bp] = group_bp(empty_struct,subj_list,setpath, setfil,task,type,ses) 
    
    subjs_bp = empty_struct;
    
    % Define the channels
    load('chanlocs.mat')
    subjs_bp(1).chanlocs = chanlocs;
    ch_names=extractfield(subjs_bp.chanlocs,'labels');
    nb_subjects = length(subj_list);
    % Define cell structure for all powers and 31 channels
    subjs_bp.PI.relative.delta = NaN(31,nb_subjects); subjs_bp.PI.relative.theta = NaN(31,nb_subjects);
    subjs_bp.PI.relative.alpha = NaN(31,nb_subjects); subjs_bp.PI.relative.beta = NaN(31,nb_subjects);

    subjs_bp.PI.absolute.delta = NaN(31,nb_subjects); subjs_bp.PI.absolute.theta = NaN(31,nb_subjects);
    subjs_bp.PI.absolute.alpha = NaN(31,nb_subjects); subjs_bp.PI.absolute.beta = NaN(31,nb_subjects);
    
    subjs_bp.PR.relative.delta = NaN(31,nb_subjects); subjs_bp.PR.relative.theta = NaN(31,nb_subjects);
    subjs_bp.PR.relative.alpha = NaN(31,nb_subjects); subjs_bp.PR.relative.beta = NaN(31,nb_subjects);

    subjs_bp.PR.absolute.delta = NaN(31,nb_subjects); subjs_bp.PR.absolute.theta = NaN(31,nb_subjects);
    subjs_bp.PR.absolute.alpha = NaN(31,nb_subjects); subjs_bp.PR.absolute.beta = NaN(31,nb_subjects);

    subjs_bp.PI.relative_sum.delta = NaN(31,nb_subjects); subjs_bp.PI.relative_sum.theta = NaN(31,nb_subjects);
    subjs_bp.PI.relative_sum.alpha = NaN(31,nb_subjects); subjs_bp.PI.relative_sum.beta = NaN(31,nb_subjects);

    subjs_bp.PI.absolute_sum.delta = NaN(31,nb_subjects); subjs_bp.PI.absolute_sum.theta = NaN(31,nb_subjects);
    subjs_bp.PI.absolute_sum.alpha = NaN(31,nb_subjects); subjs_bp.PI.absolute_sum.beta = NaN(31,nb_subjects);
    
    subjs_bp.PR.relative_sum.delta = NaN(31,nb_subjects); subjs_bp.PR.relative_sum.theta = NaN(31,nb_subjects);
    subjs_bp.PR.relative_sum.alpha = NaN(31,nb_subjects); subjs_bp.PR.relative_sum.beta = NaN(31,nb_subjects);

    subjs_bp.PR.absolute_sum.delta = NaN(31,nb_subjects); subjs_bp.PR.absolute_sum.theta = NaN(31,nb_subjects);
    subjs_bp.PR.absolute_sum.alpha = NaN(31,nb_subjects); subjs_bp.PR.absolute_sum.beta = NaN(31,nb_subjects);
    
    subjs_bp.PI.average_power.delta = NaN(31,nb_subjects); subjs_bp.PI.average_power.theta = NaN(31,nb_subjects);
    subjs_bp.PI.average_power.alpha = NaN(31,nb_subjects); subjs_bp.PI.average_power.beta = NaN(31,nb_subjects);
    
    subjs_bp.PR.average_power.delta = NaN(31,nb_subjects); subjs_bp.PR.average_power.theta = NaN(31,nb_subjects);
    subjs_bp.PR.average_power.alpha = NaN(31,nb_subjects); subjs_bp.PR.average_power.beta = NaN(31,nb_subjects);

    subjs_bp.PI.bandpower_ratio.delta_theta = NaN(31,nb_subjects);
    subjs_bp.PI.bandpower_ratio.delta_alpha = NaN(31,nb_subjects);
    subjs_bp.PI.bandpower_ratio.delta_beta = NaN(31,nb_subjects);
    subjs_bp.PI.bandpower_ratio.theta_alpha = NaN(31,nb_subjects);
    subjs_bp.PI.bandpower_ratio.theta_beta = NaN(31,nb_subjects);
    subjs_bp.PI.bandpower_ratio.alpha_beta = NaN(31,nb_subjects);

    subjs_bp.PR.bandpower_ratio.delta_theta = NaN(31,nb_subjects);
    subjs_bp.PR.bandpower_ratio.delta_alpha = NaN(31,nb_subjects);
    subjs_bp.PR.bandpower_ratio.delta_beta = NaN(31,nb_subjects);
    subjs_bp.PR.bandpower_ratio.theta_alpha = NaN(31,nb_subjects);
    subjs_bp.PR.bandpower_ratio.theta_beta = NaN(31,nb_subjects);
    subjs_bp.PR.bandpower_ratio.alpha_beta = NaN(31,nb_subjects);

    subjs_bp.PI.mean_freq = NaN(31,nb_subjects); subjs_bp.PI.rmsf = NaN(31,nb_subjects);
    subjs_bp.PR.mean_freq = NaN(31,nb_subjects); subjs_bp.PR.rmsf = NaN(31,nb_subjects);
    
    subjs_bp.PI.total_power = NaN(31,nb_subjects);
    subjs_bp.PR.total_power = NaN(31,nb_subjects);

    subjs_bp.PI.total_power_sum = NaN(31,nb_subjects);
    subjs_bp.PR.total_power_sum = NaN(31,nb_subjects);


    subjs_bp.PI.PSD = cell(nb_subjects,1);
    subjs_bp.PR.PSD = cell(nb_subjects,1);

    subjs_bp.PI.freqbins = cell(nb_subjects,1);
    subjs_bp.PR.freqbins = cell(nb_subjects,1);

    subjs_bp.chanlocs_subj = cell(nb_subjects,1);
    
    % Add the bandpowers of each subject
    for j = 1:nb_subjects
        i = subj_list(j);
        % Load data for subject i
        EEG_PI_bp = pop_loadset(strcat(sprintf(setfil,type,i,ses,task),'_PI_bp.set'),sprintf(setpath,type,i));
        EEG_PR_bp = pop_loadset(strcat(sprintf(setfil,type,i,ses,task),'_PR_bp.set'),sprintf(setpath,type,i));

        subjs_bp.PI.PSD{j} = EEG_PI_bp.psd;
        subjs_bp.PR.PSD{j} = EEG_PR_bp.psd;

        subjs_bp.PI.freqbins{j} = EEG_PI_bp.freqbins;
        subjs_bp.PR.freqbins{j} = EEG_PR_bp.freqbins;

        subjs_bp.chanlocs_subj{j} = EEG_PI_bp.chanlocs;
    
        
        for ch = 1:EEG_PI_bp.nbchan
            % Takes into account that some channels are removed and so the
            % idx is not always the same 
            idx = find(strcmp(ch_names,EEG_PI_bp.chanlocs(ch).labels));
    
            subjs_bp.PI.relative.delta(idx,j) = EEG_PI_bp.bandpower_avg.relative.delta(ch,:);
            subjs_bp.PI.relative.theta(idx,j) = EEG_PI_bp.bandpower_avg.relative.theta(ch,:);
            subjs_bp.PI.relative.alpha(idx,j) = EEG_PI_bp.bandpower_avg.relative.alpha(ch,:);
            subjs_bp.PI.relative.beta(idx,j) = EEG_PI_bp.bandpower_avg.relative.beta(ch,:);
    
            subjs_bp.PR.relative.delta(idx,j) = EEG_PR_bp.bandpower_avg.relative.delta(ch,:);
            subjs_bp.PR.relative.theta(idx,j) = EEG_PR_bp.bandpower_avg.relative.theta(ch,:);
            subjs_bp.PR.relative.alpha(idx,j) = EEG_PR_bp.bandpower_avg.relative.alpha(ch,:);
            subjs_bp.PR.relative.beta(idx,j) = EEG_PR_bp.bandpower_avg.relative.beta(ch,:);

            subjs_bp.PI.absolute.delta(idx,j) = EEG_PI_bp.bandpower_avg.absolute.delta(ch,:);
            subjs_bp.PI.absolute.theta(idx,j) = EEG_PI_bp.bandpower_avg.absolute.theta(ch,:);
            subjs_bp.PI.absolute.alpha(idx,j) = EEG_PI_bp.bandpower_avg.absolute.alpha(ch,:);
            subjs_bp.PI.absolute.beta(idx,j) = EEG_PI_bp.bandpower_avg.absolute.beta(ch,:);
    
            subjs_bp.PR.absolute.delta(idx,j) = EEG_PR_bp.bandpower_avg.absolute.delta(ch,:);
            subjs_bp.PR.absolute.theta(idx,j) = EEG_PR_bp.bandpower_avg.absolute.theta(ch,:);
            subjs_bp.PR.absolute.alpha(idx,j) = EEG_PR_bp.bandpower_avg.absolute.alpha(ch,:);
            subjs_bp.PR.absolute.beta(idx,j) = EEG_PR_bp.bandpower_avg.absolute.beta(ch,:);
            
            subjs_bp.PI.relative_sum.delta(idx,j) = EEG_PI_bp.bandpower_sum_avg.relative.delta(ch,:);
            subjs_bp.PI.relative_sum.theta(idx,j) = EEG_PI_bp.bandpower_sum_avg.relative.theta(ch,:);
            subjs_bp.PI.relative_sum.alpha(idx,j) = EEG_PI_bp.bandpower_sum_avg.relative.alpha(ch,:);
            subjs_bp.PI.relative_sum.beta(idx,j) = EEG_PI_bp.bandpower_sum_avg.relative.beta(ch,:);
    
            subjs_bp.PR.relative_sum.delta(idx,j) = EEG_PR_bp.bandpower_sum_avg.relative.delta(ch,:);
            subjs_bp.PR.relative_sum.theta(idx,j) = EEG_PR_bp.bandpower_sum_avg.relative.theta(ch,:);
            subjs_bp.PR.relative_sum.alpha(idx,j) = EEG_PR_bp.bandpower_sum_avg.relative.alpha(ch,:);
            subjs_bp.PR.relative_sum.beta(idx,j) = EEG_PR_bp.bandpower_sum_avg.relative.beta(ch,:);

            subjs_bp.PI.absolute_sum.delta(idx,j) = EEG_PI_bp.bandpower_sum_avg.absolute.delta(ch,:);
            subjs_bp.PI.absolute_sum.theta(idx,j) = EEG_PI_bp.bandpower_sum_avg.absolute.theta(ch,:);
            subjs_bp.PI.absolute_sum.alpha(idx,j) = EEG_PI_bp.bandpower_sum_avg.absolute.alpha(ch,:);
            subjs_bp.PI.absolute_sum.beta(idx,j) = EEG_PI_bp.bandpower_sum_avg.absolute.beta(ch,:);
    
            subjs_bp.PR.absolute_sum.delta(idx,j) = EEG_PR_bp.bandpower_sum_avg.absolute.delta(ch,:);
            subjs_bp.PR.absolute_sum.theta(idx,j) = EEG_PR_bp.bandpower_sum_avg.absolute.theta(ch,:);
            subjs_bp.PR.absolute_sum.alpha(idx,j) = EEG_PR_bp.bandpower_sum_avg.absolute.alpha(ch,:);
            subjs_bp.PR.absolute_sum.beta(idx,j) = EEG_PR_bp.bandpower_sum_avg.absolute.beta(ch,:);

            subjs_bp.PI.average_power.delta(idx,j) = EEG_PI_bp.average_power_avg.delta(ch,:);
            subjs_bp.PI.average_power.theta(idx,j) = EEG_PI_bp.average_power_avg.theta(ch,:);
            subjs_bp.PI.average_power.alpha(idx,j) = EEG_PI_bp.average_power_avg.alpha(ch,:);
            subjs_bp.PI.average_power.beta(idx,j) = EEG_PI_bp.average_power_avg.beta(ch,:);
            
            subjs_bp.PR.average_power.delta(idx,j) = EEG_PR_bp.average_power_avg.delta(ch,:);
            subjs_bp.PR.average_power.theta(idx,j) = EEG_PR_bp.average_power_avg.theta(ch,:);
            subjs_bp.PR.average_power.alpha(idx,j) = EEG_PR_bp.average_power_avg.alpha(ch,:);
            subjs_bp.PR.average_power.beta(idx,j) = EEG_PR_bp.average_power_avg.beta(ch,:);

            subjs_bp.PI.bandpower_ratio.delta_theta(idx,j) = EEG_PI_bp.bandpower_ratio_avg.delta_theta(ch,:);
            subjs_bp.PI.bandpower_ratio.delta_alpha(idx,j) = EEG_PI_bp.bandpower_ratio_avg.delta_alpha(ch,:);
            subjs_bp.PI.bandpower_ratio.delta_beta(idx,j) = EEG_PI_bp.bandpower_ratio_avg.delta_beta(ch,:);
            subjs_bp.PI.bandpower_ratio.theta_alpha(idx,j) = EEG_PI_bp.bandpower_ratio_avg.theta_alpha(ch,:);
            subjs_bp.PI.bandpower_ratio.theta_beta(idx,j) = EEG_PI_bp.bandpower_ratio_avg.theta_beta(ch,:);
            subjs_bp.PI.bandpower_ratio.alpha_beta(idx,j) = EEG_PI_bp.bandpower_ratio_avg.alpha_beta(ch,:);
            
            subjs_bp.PR.bandpower_ratio.delta_theta(idx,j) = EEG_PR_bp.bandpower_ratio_avg.delta_theta(ch,:);
            subjs_bp.PR.bandpower_ratio.delta_alpha(idx,j) = EEG_PR_bp.bandpower_ratio_avg.delta_alpha(ch,:);
            subjs_bp.PR.bandpower_ratio.delta_beta(idx,j) = EEG_PR_bp.bandpower_ratio_avg.delta_beta(ch,:);
            subjs_bp.PR.bandpower_ratio.theta_alpha(idx,j) = EEG_PR_bp.bandpower_ratio_avg.theta_alpha(ch,:);
            subjs_bp.PR.bandpower_ratio.theta_beta(idx,j) = EEG_PR_bp.bandpower_ratio_avg.theta_beta(ch,:);
            subjs_bp.PR.bandpower_ratio.alpha_beta(idx,j) = EEG_PR_bp.bandpower_ratio_avg.alpha_beta(ch,:);

            subjs_bp.PI.mean_freq(idx,j) = EEG_PI_bp.mean_freq_avg(ch,1);
            subjs_bp.PI.rmsf(idx,j) = EEG_PI_bp.rmsf_avg(ch,1);
            subjs_bp.PR.mean_freq(idx,j) = EEG_PR_bp.mean_freq_avg(ch,1);
            subjs_bp.PR.rmsf(idx,j) = EEG_PR_bp.rmsf_avg(ch,1);
            
            subjs_bp.PI.total_power(idx,j) = EEG_PI_bp.total_power_avg(ch,1);
            subjs_bp.PR.total_power(idx,j) = EEG_PR_bp.total_power_avg(ch,1);

            subjs_bp.PI.total_power_sum(idx,j) = EEG_PI_bp.total_power_sum_avg(ch,1);
            subjs_bp.PR.total_power_sum(idx,j) = EEG_PR_bp.total_power_sum_avg(ch,1);
        end
    
    end
    

    subjs_bp.PI_avg.relative.delta = mean(subjs_bp.PI.relative.delta,2,"omitnan");
    subjs_bp.PI_avg.relative.theta = mean(subjs_bp.PI.relative.theta,2,"omitnan");
    subjs_bp.PI_avg.relative.alpha = mean(subjs_bp.PI.relative.alpha,2,"omitnan");
    subjs_bp.PI_avg.relative.beta = mean(subjs_bp.PI.relative.beta,2,"omitnan");

    subjs_bp.PR_avg.relative.delta = mean(subjs_bp.PR.relative.delta,2,"omitnan");
    subjs_bp.PR_avg.relative.theta = mean(subjs_bp.PR.relative.theta,2,"omitnan");
    subjs_bp.PR_avg.relative.alpha = mean(subjs_bp.PR.relative.alpha,2,"omitnan");
    subjs_bp.PR_avg.relative.beta = mean(subjs_bp.PR.relative.beta,2,"omitnan");

    subjs_bp.PI_avg.absolute.delta = mean(subjs_bp.PI.absolute.delta,2,"omitnan");
    subjs_bp.PI_avg.absolute.theta = mean(subjs_bp.PI.absolute.theta,2,"omitnan");
    subjs_bp.PI_avg.absolute.alpha = mean(subjs_bp.PI.absolute.alpha,2,"omitnan");
    subjs_bp.PI_avg.absolute.beta = mean(subjs_bp.PI.absolute.beta,2,"omitnan");

    subjs_bp.PR_avg.absolute.delta = mean(subjs_bp.PR.absolute.delta,2,"omitnan");
    subjs_bp.PR_avg.absolute.theta = mean(subjs_bp.PR.absolute.theta,2,"omitnan");
    subjs_bp.PR_avg.absolute.alpha = mean(subjs_bp.PR.absolute.alpha,2,"omitnan");
    subjs_bp.PR_avg.absolute.beta = mean(subjs_bp.PR.absolute.beta,2,"omitnan");

    subjs_bp.PI_avg.relative_sum.delta = mean(subjs_bp.PI.relative_sum.delta,2,"omitnan");
    subjs_bp.PI_avg.relative_sum.theta = mean(subjs_bp.PI.relative_sum.theta,2,"omitnan");
    subjs_bp.PI_avg.relative_sum.alpha = mean(subjs_bp.PI.relative_sum.alpha,2,"omitnan");
    subjs_bp.PI_avg.relative_sum.beta = mean(subjs_bp.PI.relative_sum.beta,2,"omitnan");

    subjs_bp.PR_avg.relative_sum.delta = mean(subjs_bp.PR.relative_sum.delta,2,"omitnan");
    subjs_bp.PR_avg.relative_sum.theta = mean(subjs_bp.PR.relative_sum.theta,2,"omitnan");
    subjs_bp.PR_avg.relative_sum.alpha = mean(subjs_bp.PR.relative_sum.alpha,2,"omitnan");
    subjs_bp.PR_avg.relative_sum.beta = mean(subjs_bp.PR.relative_sum.beta,2,"omitnan");

    subjs_bp.PI_avg.absolute_sum.delta = mean(subjs_bp.PI.absolute_sum.delta,2,"omitnan");
    subjs_bp.PI_avg.absolute_sum.theta = mean(subjs_bp.PI.absolute_sum.theta,2,"omitnan");
    subjs_bp.PI_avg.absolute_sum.alpha = mean(subjs_bp.PI.absolute_sum.alpha,2,"omitnan");
    subjs_bp.PI_avg.absolute_sum.beta = mean(subjs_bp.PI.absolute_sum.beta,2,"omitnan");

    subjs_bp.PR_avg.absolute_sum.delta = mean(subjs_bp.PR.absolute_sum.delta,2,"omitnan");
    subjs_bp.PR_avg.absolute_sum.theta = mean(subjs_bp.PR.absolute_sum.theta,2,"omitnan");
    subjs_bp.PR_avg.absolute_sum.alpha = mean(subjs_bp.PR.absolute_sum.alpha,2,"omitnan");
    subjs_bp.PR_avg.absolute_sum.beta = mean(subjs_bp.PR.absolute_sum.beta,2,"omitnan");
    
    subjs_bp.PI_avg.average_power.delta = mean(subjs_bp.PI.average_power.delta,2,"omitnan");
    subjs_bp.PI_avg.average_power.theta = mean(subjs_bp.PI.average_power.theta,2,"omitnan");
    subjs_bp.PI_avg.average_power.alpha = mean(subjs_bp.PI.average_power.alpha,2,"omitnan");
    subjs_bp.PI_avg.average_power.beta = mean(subjs_bp.PI.average_power.beta,2,"omitnan");

    subjs_bp.PR_avg.average_power.delta = mean(subjs_bp.PR.average_power.delta,2,"omitnan");
    subjs_bp.PR_avg.average_power.theta = mean(subjs_bp.PR.average_power.theta,2,"omitnan");
    subjs_bp.PR_avg.average_power.alpha = mean(subjs_bp.PR.average_power.alpha,2,"omitnan");
    subjs_bp.PR_avg.average_power.beta = mean(subjs_bp.PR.average_power.beta,2,"omitnan");

    subjs_bp.PI_avg.bandpower_ratio.delta_theta = mean(subjs_bp.PI.bandpower_ratio.delta_theta,2,"omitnan");
    subjs_bp.PI_avg.bandpower_ratio.delta_alpha = mean(subjs_bp.PI.bandpower_ratio.delta_alpha,2,"omitnan");
    subjs_bp.PI_avg.bandpower_ratio.delta_beta = mean(subjs_bp.PI.bandpower_ratio.delta_beta,2,"omitnan");
    subjs_bp.PI_avg.bandpower_ratio.theta_alpha = mean(subjs_bp.PI.bandpower_ratio.theta_alpha,2,"omitnan");
    subjs_bp.PI_avg.bandpower_ratio.theta_beta = mean(subjs_bp.PI.bandpower_ratio.theta_beta,2,"omitnan");
    subjs_bp.PI_avg.bandpower_ratio.alpha_beta = mean(subjs_bp.PI.bandpower_ratio.alpha_beta,2,"omitnan");

    subjs_bp.PR_avg.bandpower_ratio.delta_theta = mean(subjs_bp.PR.bandpower_ratio.delta_theta,2,"omitnan");
    subjs_bp.PR_avg.bandpower_ratio.delta_alpha = mean(subjs_bp.PR.bandpower_ratio.delta_alpha,2,"omitnan");
    subjs_bp.PR_avg.bandpower_ratio.delta_beta = mean(subjs_bp.PR.bandpower_ratio.delta_beta,2,"omitnan");
    subjs_bp.PR_avg.bandpower_ratio.theta_alpha = mean(subjs_bp.PR.bandpower_ratio.theta_alpha,2,"omitnan");
    subjs_bp.PR_avg.bandpower_ratio.theta_beta = mean(subjs_bp.PR.bandpower_ratio.theta_beta,2,"omitnan");
    subjs_bp.PR_avg.bandpower_ratio.alpha_beta = mean(subjs_bp.PR.bandpower_ratio.alpha_beta,2,"omitnan");

    subjs_bp.PI_avg.mean_freq = mean(subjs_bp.PI.mean_freq,2,"omitnan");
    subjs_bp.PI_avg.rmsf =  mean(subjs_bp.PI.rmsf,2,"omitnan");
    subjs_bp.PR_avg.mean_freq =  mean(subjs_bp.PR.mean_freq,2,"omitnan");
    subjs_bp.PR_avg.rmsf = mean(subjs_bp.PR.rmsf,2,"omitnan");

    subjs_bp.PI_avg.total_power =  mean(subjs_bp.PI.total_power,2,"omitnan");
    subjs_bp.PR_avg.total_power = mean(subjs_bp.PR.total_power,2,"omitnan");

    subjs_bp.PI_avg.total_power_sum =  mean(subjs_bp.PI.total_power_sum,2,"omitnan");
    subjs_bp.PR_avg.total_power_sum = mean(subjs_bp.PR.total_power_sum,2,"omitnan");



      

end

%             subjs_bp.PI.delta{idx}(end+1:end+epochs) = EEG_PI_bp.bandpower.delta(j,:);
%             subjs_bp.PI.theta{idx}(end+1:end+epochs) = EEG_PI_bp.bandpower.theta(j,:);
%             subjs_bp.PI.alpha{idx}(end+1:end+epochs) = EEG_PI_bp.bandpower.alpha(j,:);
%             subjs_bp.PI.beta{idx}(end+1:end+epochs) = EEG_PI_bp.bandpower.beta(j,:);
%     
%             subjs_bp.PR.delta{idx}(end+1:end+epochs) = EEG_PR_bp.bandpower.delta(j,:);
%             subjs_bp.PR.theta{idx}(end+1:end+epochs) = EEG_PR_bp.bandpower.theta(j,:);
%             subjs_bp.PR.alpha{idx}(end+1:end+epochs) = EEG_PR_bp.bandpower.alpha(j,:);
%             subjs_bp.PR.beta{idx}(end+1:end+epochs) = EEG_PR_bp.bandpower.beta(j,:);
function EEG_bp = compute_bandpower(EEG,method,type)

% method: 'all' -> uses all epochs; 'auto' -> uses only epoch that don't have bad segments
%         'manual' -> uses only epochs that weren't marked manually as bad
% type: PI or PR

    EEG_bp = EEG;
    fs = EEG_bp.srate;
    
    switch method
        case 'all'
            for i = 1:EEG_bp.trials
                epoch_data = EEG_bp.data(:,:,i);
                epoch_data = epoch_data-mean(epoch_data,2);
            
                % Compute PSD Welch estimate: win = 2seg; overlap = 50%.
                win = 2*fs;
                [pxx,f] = pwelch(epoch_data',win,0.5*win,fs*2,fs);
                EEG_bp.freqbins = f;
                EEG_bp.psd(:,:,i) = pxx';
            
                % Compute bandpower integral for delta [1-4]; theta [4-8]; alpha [8-13]; beta [13-30] Hz
                total_power = bandpower(EEG_bp.psd(:,:,i)',f,[1,30],'psd')';
                EEG_bp.total_power(:,i) = total_power;
                EEG_bp.bandpower.relative.delta(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[1,4],'psd')'./total_power;
                EEG_bp.bandpower.relative.theta(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[4,8],'psd')'./total_power;
                EEG_bp.bandpower.relative.alpha(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[8,13],'psd')'./total_power;
                EEG_bp.bandpower.relative.beta(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[13,30],'psd')'./total_power;

                EEG_bp.bandpower.absolute.delta(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[1,4],'psd')';
                EEG_bp.bandpower.absolute.theta(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[4,8],'psd')';
                EEG_bp.bandpower.absolute.alpha(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[8,13],'psd')';
                EEG_bp.bandpower.absolute.beta(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[13,30],'psd')';
                
                % Bandpower with sum
                total_power_sum = sum(EEG_bp.psd(:,3:61,i),2);
                EEG_bp.total_power_sum(:,i) = total_power_sum;
                EEG_bp.bandpower_sum.relative.delta(:,i) = sum(EEG_bp.psd(:,3:9,i),2)./total_power_sum;
                EEG_bp.bandpower_sum.relative.theta(:,i) = sum(EEG_bp.psd(:,9:17,i),2)./total_power_sum;
                EEG_bp.bandpower_sum.relative.alpha(:,i) = sum(EEG_bp.psd(:,17:27,i),2)./total_power_sum;
                EEG_bp.bandpower_sum.relative.beta(:,i) = sum(EEG_bp.psd(:,27:61,i),2)./total_power_sum;

                EEG_bp.bandpower_sum.absolute.delta(:,i) = sum(EEG_bp.psd(:,3:9,i),2);
                EEG_bp.bandpower_sum.absolute.theta(:,i) = sum(EEG_bp.psd(:,9:17,i),2);
                EEG_bp.bandpower_sum.absolute.alpha(:,i) = sum(EEG_bp.psd(:,17:27,i),2);
                EEG_bp.bandpower_sum.absolute.beta(:,i) = sum(EEG_bp.psd(:,27:61,i),2);

                % Average Power
                EEG_bp.average_power.delta(:,i) = sum(EEG_bp.psd(:,3:9,i),2)/7;
                EEG_bp.average_power.theta(:,i) = sum(EEG_bp.psd(:,9:17,i),2)/9;
                EEG_bp.average_power.alpha(:,i) = sum(EEG_bp.psd(:,17:27,i),2)/11;
                EEG_bp.average_power.beta(:,i) = sum(EEG_bp.psd(:,27:61,i),2)/35;
                
                % Bandpower ratios
                EEG_bp.bandpower_ratio.delta_theta = EEG_bp.bandpower.relative.delta./EEG_bp.bandpower.relative.theta;
                EEG_bp.bandpower_ratio.delta_alpha = EEG_bp.bandpower.relative.delta./EEG_bp.bandpower.relative.alpha;
                EEG_bp.bandpower_ratio.delta_beta = EEG_bp.bandpower.relative.delta./EEG_bp.bandpower.relative.beta;
                EEG_bp.bandpower_ratio.theta_alpha = EEG_bp.bandpower.relative.theta./EEG_bp.bandpower.relative.alpha;
                EEG_bp.bandpower_ratio.theta_beta = EEG_bp.bandpower.relative.theta./EEG_bp.bandpower.relative.beta;
                EEG_bp.bandpower_ratio.alpha_beta = EEG_bp.bandpower.relative.alpha./EEG_bp.bandpower.relative.beta;


                % Compute Mean Frequency & Root Mean Squared Frequency (RMSF)
                [~,min_idx] = min(abs(1-EEG_bp.freqbins));
                [~,max_idx] = min(abs(30-EEG_bp.freqbins));
                EEG_bp.mean_freq(:,i) = meanfreq(EEG_bp.psd(:,min_idx:max_idx,i)',EEG_bp.freqbins(min_idx:max_idx))';
                EEG_bp.rmsf(:,i) = rmsf(EEG_bp.psd(:,min_idx:max_idx,i)',EEG_bp.freqbins(min_idx:max_idx));

            end
            % Average over all trials
            EEG_bp.total_power_avg = mean(EEG_bp.total_power,2);
            EEG_bp.bandpower_avg.relative.delta = mean(EEG_bp.bandpower.relative.delta,2);
            EEG_bp.bandpower_avg.relative.theta = mean(EEG_bp.bandpower.relative.theta,2);
            EEG_bp.bandpower_avg.relative.alpha = mean(EEG_bp.bandpower.relative.alpha,2);
            EEG_bp.bandpower_avg.relative.beta = mean(EEG_bp.bandpower.relative.beta,2);

            EEG_bp.bandpower_avg.absolute.delta = mean(EEG_bp.bandpower.absolute.delta,2);
            EEG_bp.bandpower_avg.absolute.theta = mean(EEG_bp.bandpower.absolute.theta,2);
            EEG_bp.bandpower_avg.absolute.alpha = mean(EEG_bp.bandpower.absolute.alpha,2);
            EEG_bp.bandpower_avg.absolute.beta = mean(EEG_bp.bandpower.absolute.beta,2);
            
            EEG_bp.total_power_sum_avg = mean(EEG_bp.total_power_sum,2);
            EEG_bp.bandpower_sum_avg.relative.delta = mean(EEG_bp.bandpower_sum.relative.delta,2);
            EEG_bp.bandpower_sum_avg.relative.theta = mean(EEG_bp.bandpower_sum.relative.theta,2);
            EEG_bp.bandpower_sum_avg.relative.alpha = mean(EEG_bp.bandpower_sum.relative.alpha,2);
            EEG_bp.bandpower_sum_avg.relative.beta = mean(EEG_bp.bandpower_sum.relative.beta,2);

            EEG_bp.bandpower_sum_avg.absolute.delta = mean(EEG_bp.bandpower_sum.absolute.delta,2);
            EEG_bp.bandpower_sum_avg.absolute.theta = mean(EEG_bp.bandpower_sum.absolute.theta,2);
            EEG_bp.bandpower_sum_avg.absolute.alpha = mean(EEG_bp.bandpower_sum.absolute.alpha,2);
            EEG_bp.bandpower_sum_avg.absolute.beta = mean(EEG_bp.bandpower_sum.absolute.beta,2);

            EEG_bp.average_power_avg.delta = mean(EEG_bp.average_power.delta,2);
            EEG_bp.average_power_avg.theta = mean(EEG_bp.average_power.theta,2);
            EEG_bp.average_power_avg.alpha = mean(EEG_bp.average_power.alpha,2);
            EEG_bp.average_power_avg.beta = mean(EEG_bp.average_power.beta,2);

            EEG_bp.bandpower_ratio_avg.delta_theta = mean(EEG_bp.bandpower_ratio.delta_theta,2);
            EEG_bp.bandpower_ratio_avg.delta_alpha = mean(EEG_bp.bandpower_ratio.delta_alpha,2);
            EEG_bp.bandpower_ratio_avg.delta_beta = mean(EEG_bp.bandpower_ratio.delta_beta,2);
            EEG_bp.bandpower_ratio_avg.theta_alpha = mean(EEG_bp.bandpower_ratio.theta_alpha,2);
            EEG_bp.bandpower_ratio_avg.theta_beta = mean(EEG_bp.bandpower_ratio.theta_beta,2);
            EEG_bp.bandpower_ratio_avg.alpha_beta = mean(EEG_bp.bandpower_ratio.alpha_beta,2);

            EEG_bp.mean_freq_avg =  mean(EEG_bp.mean_freq,2); 
            EEG_bp.rmsf_avg = mean(EEG_bp.rmsf,2);


        case 'manual'
            for i = 1:EEG_bp.trials
                if  EEG_bp.task.(strcat(type,'_badsamples')){2,i} == 0

                    epoch_data = EEG_bp.data(:,:,i);
                    epoch_data = epoch_data-mean(epoch_data,2);

                    % Compute PSD Welch estimate: win = 2seg; overlap = 50%.
                    win = 2*fs;
                    [pxx,f] = pwelch(epoch_data',win,0.5*win,fs*2,fs);
                    EEG_bp.freqbins = f;
                    EEG_bp.psd(:,:,i) = pxx';

                    % Compute bandpower for delta [1-4]; theta [4-8]; alpha [8-13]; beta
                    % [13-30] Hz
                    total_power = bandpower(EEG_bp.psd(:,:,i)',f,[1,30],'psd')';
                    EEG_bp.total_power(:,i) = total_power;
                    EEG_bp.bandpower.relative.delta(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[1,4],'psd')'./total_power;
                    EEG_bp.bandpower.relative.theta(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[4,8],'psd')'./total_power;
                    EEG_bp.bandpower.relative.alpha(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[8,13],'psd')'./total_power;
                    EEG_bp.bandpower.relative.beta(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[13,30],'psd')'./total_power;

                    EEG_bp.bandpower.absolute.delta(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[1,4],'psd')';
                    EEG_bp.bandpower.absolute.theta(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[4,8],'psd')';
                    EEG_bp.bandpower.absolute.alpha(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[8,13],'psd')';
                    EEG_bp.bandpower.absolute.beta(:,i) = bandpower(EEG_bp.psd(:,:,i)',f,[13,30],'psd')';
                    
                    % Bandpower with sum
                    total_power_sum = sum(EEG_bp.psd(:,3:61,i),2);
                    EEG_bp.total_power_sum(:,i) = total_power_sum;
                    EEG_bp.bandpower_sum.relative.delta(:,i) = sum(EEG_bp.psd(:,3:9,i),2)./total_power_sum;
                    EEG_bp.bandpower_sum.relative.theta(:,i) = sum(EEG_bp.psd(:,9:17,i),2)./total_power_sum;
                    EEG_bp.bandpower_sum.relative.alpha(:,i) = sum(EEG_bp.psd(:,17:27,i),2)./total_power_sum;
                    EEG_bp.bandpower_sum.relative.beta(:,i) = sum(EEG_bp.psd(:,27:61,i),2)./total_power_sum;

                    EEG_bp.bandpower_sum.absolute.delta(:,i) = sum(EEG_bp.psd(:,3:9,i),2);
                    EEG_bp.bandpower_sum.absolute.theta(:,i) = sum(EEG_bp.psd(:,9:17,i),2);
                    EEG_bp.bandpower_sum.absolute.alpha(:,i) = sum(EEG_bp.psd(:,17:27,i),2);
                    EEG_bp.bandpower_sum.absolute.beta(:,i) = sum(EEG_bp.psd(:,27:61,i),2);
                    % 

                    EEG_bp.average_power.delta(:,i) = sum(EEG_bp.psd(:,3:9,i),2)/7;
                    EEG_bp.average_power.theta(:,i) = sum(EEG_bp.psd(:,9:17,i),2)/9;
                    EEG_bp.average_power.alpha(:,i) = sum(EEG_bp.psd(:,17:27,i),2)/11;
                    EEG_bp.average_power.beta(:,i) = sum(EEG_bp.psd(:,27:61,i),2)/35;

                    % Compute Mean Frequency & Root Mean Squared Frequency (RMSF)
                    [~,min_idx] = min(abs(1-EEG_bp.freqbins));
                    [~,max_idx] = min(abs(30-EEG_bp.freqbins));
                    EEG_bp.mean_freq(:,i) = meanfreq(EEG_bp.psd(:,min_idx:max_idx,i)',EEG_bp.freqbins(min_idx:max_idx))';
                    EEG_bp.rmsf(:,i) = rmsf(EEG_bp.psd(:,min_idx:max_idx,i)',EEG_bp.freqbins(min_idx:max_idx));
                    
                else
                    EEG_bp.psd(:,:,i) = NaN(EEG_bp.nbchan,251);

                    EEG_bp.total_power(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower.relative.delta(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower.relative.theta(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower.relative.alpha(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower.relative.beta(:,i) = NaN(EEG_bp.nbchan,1);

                    EEG_bp.bandpower.absolute.delta(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower.absolute.theta(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower.absolute.alpha(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower.absolute.beta(:,i) = NaN(EEG_bp.nbchan,1);

                    EEG_bp.total_power_sum(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower_sum.relative.delta(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower_sum.relative.theta(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower_sum.relative.alpha(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower_sum.relative.beta(:,i) = NaN(EEG_bp.nbchan,1);

                    EEG_bp.bandpower_sum.absolute.delta(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower_sum.absolute.theta(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower_sum.absolute.alpha(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.bandpower_sum.absolute.beta(:,i) = NaN(EEG_bp.nbchan,1);   
                    
                    EEG_bp.average_power.delta(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.average_power.theta(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.average_power.alpha(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.average_power.beta(:,i) = NaN(EEG_bp.nbchan,1);

                    EEG_bp.mean_freq(:,i) = NaN(EEG_bp.nbchan,1);
                    EEG_bp.rmsf(:,i) = NaN(EEG_bp.nbchan,1);

                end
            end
                % Compute bandpower ratios
                EEG_bp.bandpower_ratio.delta_theta = EEG_bp.bandpower.relative.delta./EEG_bp.bandpower.relative.theta;
                EEG_bp.bandpower_ratio.delta_alpha = EEG_bp.bandpower.relative.delta./EEG_bp.bandpower.relative.alpha;
                EEG_bp.bandpower_ratio.delta_beta = EEG_bp.bandpower.relative.delta./EEG_bp.bandpower.relative.beta;
                EEG_bp.bandpower_ratio.theta_alpha = EEG_bp.bandpower.relative.theta./EEG_bp.bandpower.relative.alpha;
                EEG_bp.bandpower_ratio.theta_beta = EEG_bp.bandpower.relative.theta./EEG_bp.bandpower.relative.beta;
                EEG_bp.bandpower_ratio.alpha_beta = EEG_bp.bandpower.relative.alpha./EEG_bp.bandpower.relative.beta;

                % Average over all trials
                EEG_bp.total_power_avg = mean(EEG_bp.total_power,2,'omitnan');
                EEG_bp.bandpower_avg.relative.delta = mean(EEG_bp.bandpower.relative.delta,2,'omitnan');
                EEG_bp.bandpower_avg.relative.theta = mean(EEG_bp.bandpower.relative.theta,2,'omitnan');
                EEG_bp.bandpower_avg.relative.alpha = mean(EEG_bp.bandpower.relative.alpha,2,'omitnan');
                EEG_bp.bandpower_avg.relative.beta = mean(EEG_bp.bandpower.relative.beta,2,'omitnan');

                EEG_bp.bandpower_avg.absolute.delta = mean(EEG_bp.bandpower.absolute.delta,2,'omitnan');
                EEG_bp.bandpower_avg.absolute.theta = mean(EEG_bp.bandpower.absolute.theta,2,'omitnan');
                EEG_bp.bandpower_avg.absolute.alpha = mean(EEG_bp.bandpower.absolute.alpha,2,'omitnan');
                EEG_bp.bandpower_avg.absolute.beta = mean(EEG_bp.bandpower.absolute.beta,2,'omitnan');
                
                EEG_bp.total_power_sum_avg = mean(EEG_bp.total_power_sum,2,'omitnan');
                EEG_bp.bandpower_sum_avg.relative.delta = mean(EEG_bp.bandpower_sum.relative.delta,2,'omitnan');
                EEG_bp.bandpower_sum_avg.relative.theta = mean(EEG_bp.bandpower_sum.relative.theta,2,'omitnan');
                EEG_bp.bandpower_sum_avg.relative.alpha = mean(EEG_bp.bandpower_sum.relative.alpha,2,'omitnan');
                EEG_bp.bandpower_sum_avg.relative.beta = mean(EEG_bp.bandpower_sum.relative.beta,2,'omitnan');

                EEG_bp.bandpower_sum_avg.absolute.delta = mean(EEG_bp.bandpower_sum.absolute.delta,2,'omitnan');
                EEG_bp.bandpower_sum_avg.absolute.theta = mean(EEG_bp.bandpower_sum.absolute.theta,2,'omitnan');
                EEG_bp.bandpower_sum_avg.absolute.alpha = mean(EEG_bp.bandpower_sum.absolute.alpha,2,'omitnan');
                EEG_bp.bandpower_sum_avg.absolute.beta = mean(EEG_bp.bandpower_sum.absolute.beta,2,'omitnan');

                EEG_bp.average_power_avg.delta = mean(EEG_bp.average_power.delta,2,'omitnan');
                EEG_bp.average_power_avg.theta = mean(EEG_bp.average_power.theta,2,'omitnan');
                EEG_bp.average_power_avg.alpha = mean(EEG_bp.average_power.alpha,2,'omitnan');
                EEG_bp.average_power_avg.beta = mean(EEG_bp.average_power.beta,2,'omitnan');

                EEG_bp.bandpower_ratio_avg.delta_theta = mean(EEG_bp.bandpower_ratio.delta_theta,2,'omitnan');
                EEG_bp.bandpower_ratio_avg.delta_alpha = mean(EEG_bp.bandpower_ratio.delta_alpha,2,'omitnan');
                EEG_bp.bandpower_ratio_avg.delta_beta = mean(EEG_bp.bandpower_ratio.delta_beta,2,'omitnan');
                EEG_bp.bandpower_ratio_avg.theta_alpha = mean(EEG_bp.bandpower_ratio.theta_alpha,2,'omitnan');
                EEG_bp.bandpower_ratio_avg.theta_beta = mean(EEG_bp.bandpower_ratio.theta_beta,2,'omitnan');
                EEG_bp.bandpower_ratio_avg.alpha_beta = mean(EEG_bp.bandpower_ratio.alpha_beta,2,'omitnan');

                EEG_bp.mean_freq_avg =  mean(EEG_bp.mean_freq,2,'omitnan');
                EEG_bp.rmsf_avg = mean(EEG_bp.rmsf,2,'omitnan');
           
        otherwise
            disp('Not a valid method')

    end

end

%plot(f,10*log10(EEG_PSD.psd(1,:,1)),LineWidth=1.5)

function [psd_rmsf] = rmsf(psd,freq)

    total_power = sum(psd);
    normalized_psd = psd./total_power;
    psd_rmsf = sqrt(sum(normalized_psd.*(freq.^2)))';

end
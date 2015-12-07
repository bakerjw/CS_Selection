function [ sumDn ] = KS_stat( periods, emp_cdf, sampleSmall, means, stdevs )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Dn = zeros(length(periods),1);
for h = 1:length(periods)
    % Sort the lnSa values at each period and calculate the
    % normal CDF
    sortedlnSa = [min(sampleSmall(:,h)); sort(sampleSmall(:,h))];
    norm_cdf = normcdf(sortedlnSa,means(h),stdevs(h));
    
    % Calculate the Dn value
    Dn(h) = max(abs(emp_cdf'-norm_cdf));
    sumDn = sum(Dn);
end

end
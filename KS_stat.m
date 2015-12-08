function [ sumDn ] = KS_stat( periods, emp_cdf, sampleSmall, means, stdevs )
% calculate sum of all KS-test statistics 

sortedlnSa = [min(sampleSmall); sort(sampleSmall)];
norm_cdf = normcdf(sortedlnSa,repmat(means,size(sampleSmall,1)+1,1),repmat(stdevs,size(sampleSmall,1)+1,1));
Dn = max(abs(repmat(emp_cdf',1,length(periods)) - norm_cdf));
sumDn = sum(Dn);

end
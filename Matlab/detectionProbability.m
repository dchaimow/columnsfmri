function p = detectionProbability(cnr,N)
% DETECTIONPROBABILITY calculates probability to detect a response.
%   p = DETECTIONPROBABILITY(cnr,N) calculated the probability to detect a
%   uni- or multivariate response in N voxels if expected contrast-to-noise
%   ratio is cnr.

if numel(cnr)>1
    p = zeros(size(cnr));
    if numel(N)==1
        N = ones(size(cnr)) * N;
    end
    for z=1:numel(cnr)
        p(z) = detectionProbability(cnr(z),N(z));
    end
else    
    alpha = 0.05;
    x_crit = chi2inv(1-alpha,N);
    a = N/2;
    b = 2*(1+cnr^2);
    p = gamcdf(x_crit,a,b,'upper');
end
% gammacdf results from the following prior integration:
% p = integral2(@f,0,Inf,x_crit,Inf)/(cnr^2);
%
%     function y = f(lambda,x)
%         y = ncx2pdf(x,N,lambda).*chi2pdf(lambda./cnr.^2,N);
%     end
end


function sigma = noiseModel(V,noiseType,TR,nT,differentialFlag)
% noiseModel predicts the fMRI noise level.
%   sigma = noiseModel(V,noiseType,TR,nT,differentialFlag) calculates the
%   standard deviation of fMRI noise relative to signal after differential 
%   analysis of multiple measurements nT (see below) using voxel volume V 
%   and TR.
%
%   nT is the number of volumes to be averaged:
%   it is used to scale the thermal noise factor, assuming thermal noise is
%   uncorrelated in time
%   AND it is used to scale the physiological noise factor under the
%   assuption that physiological noise is a AR(1) process with
%   q = exp(-TR/tau), tau = 15 s (Purdon and Weisskoff 1998)
%
%   with set differential flag nT/2 volumes belong to condition A and nT/2
%   volumes to condition B
%
%   The noise model is based on Triantafyllou et al. 2005.
%   It is specified as  
%     noiseType = {'3T','7T','Thermal','Physiological'}
%    or using model parameters:
%     noisetype = [k,l,T1] corresponding to kappa and lambda in 
%   Tiantafyllou et al. 2005 and time constant T1

TR0 = 5.4;
if isnumeric(noiseType)
    kappa = noiseType(1);
    lambda = noiseType(2);
    T1 = noiseType(3);
elseif strcmp(noiseType,'3T') % from estimateNoiseModelFromTriantafyllou2005.m
    kappa = 6.6567;
    lambda = 0.0129;
    T1 = 1.607;
elseif strcmp(noiseType,'7T')
    kappa = 9.9632;
    lambda = 0.0113;
    T1 = 1.939;
elseif strcmp(noiseType,'Thermal')
    kappa = 9.9632;
    lambda = 0;
    T1 = 1.939;
elseif strcmp(noiseType,'Physiological')
    kappa = 1/eps;
    lambda = 0.0113;
    T1 = 1.939;
end

if ~differentialFlag && nT~=1
    error('for multiple measurements only differential implemented!');
elseif nT==1
    F = sqrt(tanh(TR/(2*T1))/tanh(TR0/(2*T1)));
    kappa = kappa * F;   
    sigma = sqrt(1+lambda^2*kappa^2*V.^2)./(kappa.*V);
    return;
end

s = 0;
for t1 = 1:nT/2
    for t2 = 1:nT/2
        s = s + exp((-TR*abs(t1-t2))/15);
    end
end

F = sqrt(tanh(TR/(2*T1))/tanh(TR0/(2*T1)));
kappa = kappa * F;
sigma = sqrt((4./(kappa^2*V.^2*nT)) + ((2*lambda^2)/(nT/2)^2) * s);
end
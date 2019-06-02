function [neuronal,prefMap, FORIENT] = ...
sim_columnPattern(sim,rho,delta,gwnoise)
% SIM_COLUMNPATTERN simulates realistic pattern of cortical columns. 
%   [neuronal,FORIENT] = SIM_COLUMNPATTERN(sim,rho,delta,gwnoise) simulates
%   the differential neuronal response of a pattern of cortical columns by
%   filtering of spatial Gaussian white noise gwnoise using a spatial
%   band-pass filter parameterized by main spatial frequency rho and
%   relative irregularity delta. The simulation grid is defined in sim. 
%   The simulated pattern is returned in neuronal and the spatial 
%   representation of the filter in FORIENT.
%
%   see also setupsim, sim_gwnoise

fwhmfactor = 2*sqrt(2*log(2));
sim.r = sqrt(sim.k1.^2+sim.k2.^2);
if delta == 0
    FORIENTNotNormalized = abs(sim.r - rho)<sim.dk/2;
else
    % sum of Gaussians (positive and negative):
    FORIENTNotNormalized = normpdf(sim.r,rho,(rho * delta)/fwhmfactor) + ...
        normpdf(sim.r,-rho,(rho * delta)/fwhmfactor);
end
% normalization, such that single condition model produces average response
% of 1
C = (sqrt(meanpower(FORIENTNotNormalized)))*sqrt(pi/8);
FORIENT = FORIENTNotNormalized/C;
noiseF = sim_ft2(gwnoise,sim.dx);
gamma = sim_ift2(FORIENT.*noiseF,sim.dk);
neuronal = real(gamma);
prefMap = angle(gamma)/2;




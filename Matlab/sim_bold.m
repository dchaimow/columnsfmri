function [by,psf,MTF]=sim_bold(sim,fwhm,beta,y)
% SIM_BOLD simulates spatial BOLD response. 
%   [by,PSF,MTF] = SIM_BOLD(sim,fwhm,beta,y) simulates spatial BOLD
%   response to neuronal response pattern y using a BOLD PSF with
%   full-width at half-maximum fwhm, response amplitude beta and simulation
%   grid defined in sim. The resulting BOLD response pattern is returned in
%   by. The point-spread function and modulation-transfer function are 
%   returned in PSF and MTF, respectively. 
%
%   see also setupsim

if fwhm == 0
    by = beta * y;
    psf = [];
    MTF= ones(size(y));
else
    fwfactor = 2*sqrt(2*log(2));
    psf = beta * normpdf(sim.x1,0,fwhm/fwfactor).*...
        normpdf(sim.x2,0,fwhm/fwfactor);
    MTF = beta* exp(-(sim.k1.^2+sim.k2.^2)*2*(fwhm/fwfactor)^2*pi^2);
    by = sim_ift2(MTF.*sim_ft2(y,sim.dx),sim.dk);
end
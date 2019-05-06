function [my,vWidth] = sim_mri(sim,downFactor,y)
% SIM_MRI simulates MRI voxel sampling. 
%   [my,vWidth] = SIM_MRI(sim,downFactor,y) simulates MRI voxel sampling 
%   from pattern y by reducing the k-space represenation according to 
%   downFactor and using simulation grid defined in sim. The array of 
%   sample voxels is returned in my and the voxel width is returned in 
%   vWidth.
%
%   see also setupsim

L = sim.L;
% here N is defined, so that 2N covers the FOV in one dim
N = ceil(sim.N * downFactor * 0.5);
vWidth = L/(2*N);

centerk = find(sim.k == 0);

yk = sim_ft2(y,sim.dx);
downSampledY = yk(centerk-N:centerk+N-1,centerk-N:centerk+N-1);
my = real(fftshift(ifft2(ifftshift(downSampledY))))/vWidth^2;
end



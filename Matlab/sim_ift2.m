function y = sim_ift2(fy,dk)
% SIM_IFT2 simulates 2D inverse fourier transform numerically using ifft.
%   y = SIM_IFT2(fy,dk) calculates 2D inverse Fourier transform of 2D 
%   frequency space representation fy assuming a grid spacing of dk by dk.

N = size(fy,1);
y = N^2 * dk^2 * fftshift(ifft2(ifftshift(fy)));
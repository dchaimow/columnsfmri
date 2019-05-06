function fy = sim_ft2(y,dx)
% SIM_FT2 simulates 2D fourier transform numerically using fft.
%   fy = SIM_FT2(y,dx) calculates 2D Fourier transform of 2D pattern
%   represented by matrix y assuming a grid spacing of dx by dx.

N = size(y,1);
L = N * dx;
fy = (L^2/N^2)*fftshift(fft2(ifftshift(y)));
function r = patternCorrelation(sim,orig,mri)
% PATTERNCORRELATION calculates pattern correlation. 
%   r = PATTERNCORRELATION(sim,orig,mri) calculates correlation r between
%   original pattern orig and upsampled version of pattern mri. Simulation
%   grid is defined in sim. 
%
%   see also setupsim
mriUp = real(dc_upsample(mri,sim));
C = corrcoef(orig(:),mriUp(:));        
r = C(1,2);
end

function upPattern = dc_upsample(pattern,sim)
Fy = fft2(pattern);
[m,n] = size(Fy); % size of the original matrix
nzsr = sim.N - m; % number of zero rows to add
nzsc = sim.N - n; % number of zero columns to add
% quadrants of the FFT, starting in the upper left
q1 = Fy(1:m/2,1:n/2);
q2 = Fy(1:m/2,n/2+1:n);
q3 = Fy(m/2+1:m,n/2+1:n);
q4 = Fy(m/2+1:m,1:n/2);
zpdr = zeros(nzsr,n/2);  % zeropad rows to insert
zpdc = zeros(nzsr+m,nzsc); % zeropad columns to insert
zpdFy = [ [q1;zpdr;q4] , zpdc , [q2;zpdr;q3] ]; % insert the zeros
upPattern = real(ifft2(zpdFy)) * sim.N^2/(m*n) ;
end
function sim = setupsim(N,L)
% SETUPSIM sets up simulation grid. 
%   sim = SETUPSIM(N,L) returns structure sim that defines simulation grid
%   specified by parameters:
%   N   simulation grid points along one dimension
%   L   simulation grid length along one dimension [mm]
sim.N = N;
sim.L = L;
sim.dx = L/sim.N;
sim.dk = 1/sim.L;
sim.Fs = 1/sim.dx;

sim.x = fftshift([0:sim.N/2-1 -sim.N/2:-1]*sim.dx);
sim.k = fftshift([0:sim.N/2-1 -sim.N/2:-1]*sim.dk);

[sim.x1,sim.x2] = ndgrid(sim.x,sim.x);
[sim.k1,sim.k2] = ndgrid(sim.k,sim.k);

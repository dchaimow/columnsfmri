function y = sim_gwnoise(sim)
% SIM_GWNOISE generates Gaussian white noise array. 
%   y = SIM_GWNOISE(sim) generates Gaussian white noise array according to
%   simulation grid defined in sim, to be used for simulating a 
%   column pattern using sim_columnPattern.
%
%   see also sim_columnPattern, setupsim

y = randn(sim.N,sim.N) + 1j*randn(sim.N,sim.N);


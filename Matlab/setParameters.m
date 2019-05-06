function parameters = setParameters(varargin)
% setParameters sets parameters for simulatefMRIOfColumnPatterns.
%   parameters = setParameters returns default parameters.
%
%   parameters = setParameters(s1,...) sets parameters according to
%   predefined scenarios.
%   s1,... are strings that select one of multiple scanner/pulse
%   sequence and/or pattern irregularity scenarios:
%   '3TGE','7TGE','7TSE'
%   'regular','irregular;
%
%   parameters is a structure consisting of the following fields:
%
%   -- randomNumberSeed     random number seed
%
%   -- nTrials              number of simulation trials
%   -- N                    simulation grid points along one dimension
%   -- L                    simulation grid length along one dimension [mm]
%   -- downFactors          list of k-space reduction factors for MRI
%                           simulation, determines simulated voxel widths
%
%   -- rho                  main pattern frequency, corresponds to
%                           1/(2*column spacing) [1/mm]
%   -- delta                (relative) irregularity (= bandpass FWHM/rho)
%
%   -- fwhm                 BOLD PSF FWHM [mm]
%   -- beta                 BOLD PSF amplitude
%
%   -- sliceThickness       slice thickness [mm]
%   -- AFlat                FOV area, determines number of voxels [mm^2]
%   -- TR                   TR (repetition time)
%   -- nT                   number of volumes (measurements), sum of 2
%                           conditions
%   -- noiseType            noise type, either '3T' or '7T',
%                           can alternatively be a numeric vector
%                           [kappa, lambda, T1] consisting of noise model
%                           parameters kappa and lambda (see Triantafyllou
%                           et al. 2005) and longitudinal relaxation time
%                           T1
%
%   see also simulatefMRIOfColumnPatterns

% seed for random number generator
parameters.randomNumberSeed = 23;

% simulation parameters
parameters.nTrials = 32; % number of simulation trials
parameters.N = 512;% simulation grid points
parameters.L = 24; % mm, for simulation/calculation of contrast range
% range of downsampleFactos, determines simulated voxel widths
parameters.downFactors = (2/parameters.N) * ...
    [3:23 24 25 26 27 28 30 32 34 37 40 44 48 53 60 69 80 96 120 160 240];

% pattern parameters
parameters.rho = 1/1.6; % main pattern frequency = 1/(2*column spacing)
parameters.delta = 0.5; % (moderate) irregularity (bandpass FWHM/rho)

% PSF parameters
parameters.fwhm = 1.02;  % PSF FWHM in mm
parameters.beta = 0.035; % amplitute

% further MRI parameters addecting SNR calculation
parameters.sliceThickness = 2.5; % mm
parameters.AFlat = 87; % mm^2 area for calculating nVoxels (from odc rois)
parameters.TR = 2; % s
parameters.nT = 1000; % number of volumes (sum of both conditions)
parameters.noiseType = '7T'; % alternatively: '3T'

if any(strcmp('3TGE',varargin'))
    assert(...
        ~any(strcmp('7TGE',varargin'))&&~any(strcmp('7TSE',varargin')),...
        'Only one field strength/pulse sequence scenario allowed!');
    parameters.fwhm = 2.8;  % PSF FWHM in mm
    parameters.beta = 0.03; % amplitute
    parameters.noiseType = '3T'; % alternatively: '3T'
end

if any(strcmp('7TGE',varargin'))
    assert(...
        ~any(strcmp('3TGE',varargin'))&&~any(strcmp('7TSE',varargin')),...
        'Only one field strength/pulse sequence scenario allowed!');
    parameters.fwhm = 1.02;  % PSF FWHM in mm
    parameters.beta = 0.035; % amplitute
    parameters.noiseType = '7T'; % alternatively: '3T'
end

if any(strcmp('7TSE',varargin'))
    assert(...
        ~any(strcmp('3TGE',varargin'))&&~any(strcmp('7TGE',varargin')),...
        'Only one field strength/pulse sequence scenario allowed!');
    parameters.fwhm = 0.82;  % PSF FWHM in mm
    parameters.beta = 0.025; % amplitute
    parameters.noiseType = '7T'; % alternatively: '3T'
end

if any(strcmp('irregular',varargin'))
    assert(~any(strcmp('regular',varargin')),...
        'Only one regularity scenario allowed!');
    parameters.delta = 1;
end

if any(strcmp('regular',varargin'))
    assert(~any(strcmp('irregular',varargin')),...
        'Only one regularity scenario allowed!');
    parameters.delta = 0;
end
end
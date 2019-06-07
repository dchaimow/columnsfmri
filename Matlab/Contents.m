% COLUMNSFMRI
%
% This MATLAB code supplements 
% 
% Chaimow, D., Ugurbil, K., Shmuel, A., 2018. 
% Optimization of functional MRI for detection, decoding and high-resolution 
% imaging of the response patterns of cortical columns. 
% Neuroimage. doi:10.1016/j.neuroimage.2017.04.011
% 
% This code allows to simulate fMRI of cortical column patterns using 
% different pattern and measurement parameters. It quantifies measures 
% related to possible analysis objectives and approaches, and estimates
% optimal voxel widths.
% 
% For information on how to run the model, type:
% help simulatefMRIOfColumnPatterns
%
% Files
%   simulatefMRIOfColumnPatterns - simulatefMRIOfColumnPatterns simulates fMRI of column patterns. 
%   setParameters                - setParameters sets parameters for simulatefMRIOfColumnPatterns.
%   setupsim                     - sets up simulation grid. 
%   sim_gwnoise                  - generates Gaussian white noise array. 
%   sim_columnPattern            - simulates realistic pattern of cortical columns. 
%   sim_bold                     - simulates spatial BOLD response. 
%   sim_mri                      - simulates MRI voxel sampling. 
%   detectionProbability         - calculates probability to detect a response.
%   meanpower                    - calculates the mean power.
%   noiseModel                   - noiseModel predicts the fMRI noise level.
%   patternCorrelation           - calculates pattern correlation. 
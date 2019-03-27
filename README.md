# MOLTO-3BP
Multi-Objective Low-Thrust Optimizer for Three Body Problem
[![DOI](https://zenodo.org/badge/169291732.svg)](https://zenodo.org/badge/latestdoi/169291732)



MOLTO-3BP is a fully automated Matlab tool for the preliminary design of low-thrust trajectories in the restricted three body problem. The software tool applies only for transfers departing from a Earth's circular orbit to a Moon's circular Orbit.


## Installation Guide
Installation requires simply that you download [MOLTO-3BP](https://github.com/uc3m-aerospace/MOLTO-3BP) and add all the folders and subfolders to your Matlab path.

### Dependencies
A recent version of Matlab is needed to run the code (R2016a or newer). 

## Acknowledgment

This code was developed by Nereida Agüera during her Master Thesis. Andrés Marco modified and enhanced part of the code in his Master Thesis. Thanks also to Mick Wijnen for his contribution. 


## Quick Usage Guide

In order to optimize a mission, the user needs to call the main function *molto_3bp.m* providing an input structure. The user only need to define the departure and target orbit radiuses. Here you have an example:

```matlab
% MISSION from Geostationary Earth Orbit to Low Moon Orbit.
%% Initial Data
R_e  = 6378;                                   % [km] Mean Earth's radius
R_m  = 1738;                                   % [km] Mean Moon's radius
h_i  = 35786;                                  % [km] Altitude of departure orbit
h_f  = 4000-R_m;                               % [km] Altitude of arrival orb
%
r0   = h_i+R_e;                                % [km] Distance to Earth's Center
rf   = h_f+R_m;                                % [km] Distance to Moon's Center
%
% Save to input structure
%
input.r0  = r0;
input.rf  = rf;
%
% Genetic algorithm parameters
%
input.pop = 1000;     % NSGA-II population number
input.gen = 200;      % NSGA-II generation number
input.init_file     = []; % Load initial populations file (if available)
input.output_file   = 'Example1.txt'; % Output file name
input.useParallel   = 'yes'; % (yes/no) Compute fitness function in paralle
        
% RUN MOLTO-IT ALGORITHM
output = molto_3bp(input)

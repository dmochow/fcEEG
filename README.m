%
% This toolbox implements a technique for measuring functional connectivity
% in the EEG based on spatiospectral correlations
%
% The data is first decomposed into a set of spatial components - this can
% be either PCA or ICA.  The data can also be left in the electrode domain.
%
% The component time series are then bandpassed filtered to each of the
% canonical EEG frequency bands: delta (1-4 Hz), theta (4-8), alpha (8-13), and beta (13-30).  There is
% also an option to include the gamma band (30-55).  
% Bandpassing is performed on each spatial component, such that a 6
% component time series produces a new array of 6 x 4 = 24 dimensions,
% since we use four canonical EEG frequency bands
%
% In each of these (i.e., 24) dimensions, we then compute the Hilbert
% transform and extract the magnitude time series.
%
% We then measure the functional connectivity matrix (i.e., Pearson
% correlation coefficients) among all pairs of the 24 Hilbert magnitude
% time series
%
% The result is returned as a flattened FC matrix where the upper diagonal
% entries of the correlation matrix are returned.
%
% The toolbox implements dynamic functional connectivity (dFC), which computes FC
% along successive temporal windows. 
%
% Two demonstrations are provided, both using the data from the following
% paper:
% 
% Dmochowski, J. P., Bezdek, M. A., Abelson, B. P., Johnson, J. S., Schumacher, E. H., & Parra, L. C. (2014). Audience preferences are predicted by temporal reliability of neural processing. Nature communications, 5(1), 1-9.
%
% download data from:
% https://www.dropbox.com/sh/ml8bsmr5v2dq1sg/AAAbdk3tqM1SOaFKbNONYQHBa?dl=0
% put it in ../data/superbowl or change the path below
%
% demoSlidingWindow.m demonstrates how to compute dFC
% demoStaticFc.m demonstrates how to compute static FC

% (c) Jacek P. Dmochowski, 2020

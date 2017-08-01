% to compile: mcc -m analyse_power_law.m
% to run: ./analyse_power_law.sh /usr/local/MATLAB/R2014a path_to_Ca
% -- or --
% to run: matlab -r "analyse_power_law(path_to_Ca); quit"
function rc = get_power_law(path_to_Ca,nSigma)

if nargin < 2
    nSigma = 2; % number of sigma for significant threshold
end

Ca_input = dlmread(path_to_Ca);

Ca_bi = binarise_trace(Ca_input,'None');

analyse_raster(Ca_bi,nSigma);

rc = 0;

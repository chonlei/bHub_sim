% to compile: mcc -m analyse_power_law.m
% to run: ./analyse_power_law.sh /usr/local/MATLAB/R2014a path_to_Ca
% -- or --
% to run: matlab -r "analyse_power_law(path_to_Ca); quit"
function rc = get_power_law(path_to_Ca,nSigma)

if nargin < 2
    nSigma = 2; % number of sigma for significant threshold
end

try
    Ca_input = dlmread(path_to_Ca);
catch
    Ca_input = fread(fopen(path_to_Ca, 'rb'),'float64');
    [pathstr,name,ext] = fileparts(path_to_Ca);
    fileName = [name ext];
    fileIdx = regexp(fileName,'.*model_(\d+)_morphology_(\d+)_seed_(\d+)_mode_(\d+)_.*','tokens');
    fileId = sprintf('model_%s_morphology_%s_seed_%s_mode_%s', fileIdx {1}{1}, fileIdx {1}{2}, fileIdx {1}{3}, fileIdx {1}{4});
    shape = regexp(fileName,'.*_(\w+)x(\w+)\.dat','tokens');
    Ca_input = reshape(Ca_input,str2num(shape{1}{2}),str2num(shape{1}{1}))'; % oddly they are flipped
    hubList = dlmread(fullfile(pathstr,[fileId '.log']),',');
end

Ca_bi = binarise_trace(Ca_input,'None');

analyse_raster(Ca_bi,nSigma);

rc = 0;

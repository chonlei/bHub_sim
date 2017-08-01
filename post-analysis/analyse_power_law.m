function rc = analyse_power_law(path_to_Ca)

Ca_input = dlmread(path_to_Ca);

Ca_bi = binarise_trace(Ca_input,'None');

nSigma = 2; % number of sigma for significant threshold

analyse_raster(Ca_bi,nSigma);

rc = 0;
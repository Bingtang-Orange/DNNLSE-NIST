% SETUP Initializes search paths and dependencies for DNNLSE-NIST.
% Execute this script once before running any example or test script.

fprintf('Configuring DNNLSE-NIST environment...\n');

root_dir = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(root_dir, 'src')));
addpath(fullfile(root_dir, 'example'));
addpath(fullfile(root_dir, 'test_and_data'));

if ~isempty(which('chebfun'))
    fprintf('>> Chebfun dependency detected successfully.\n');
else
    warning('Chebfun was not found on the MATLAB path. Please install Chebfun from https://www.chebfun.org/.');
end

fprintf('>> All subdirectories added to MATLAB search path. Environment ready.\n');

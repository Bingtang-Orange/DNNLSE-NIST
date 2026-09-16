% SETUP Initializes search paths and dependencies for the DNNLSE-NIST package.
% Execute this script once before running any test or example driver script.

fprintf('Configuring DNNLSE-NIST environment...\n');

root_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root_dir, 'src')));
addpath(fullfile(root_dir, 'examples'));
addpath(fullfile(root_dir, 'tests'));

% Verify Chebfun installation
if exist('chebfun', 'file') == 2
    fprintf('>> Chebfun dependency detected successfully.\n');
else
    warning('Chebfun was not found on the MATLAB path. Please install Chebfun from https://www.chebfun.org/.');
end

fprintf('>> All subdirectories added to search path. Ready to execute.\n');
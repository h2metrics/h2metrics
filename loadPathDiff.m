%% loadPathDiff
% Loads a path and associated splineData from a file.
%
% Input
%   filename
%       File to load
%
% Optional parameters
%   'autocomplete' = {true, false}
%       Finds files of the form 'filename*.mat'
%   'workdir'
%       Appends an optional working directory in front of filename
%
% Output
%   dPath
%       The path
%   phiPath
%       Path of diffeomorphisms
%   splineData (optional)
%       splineData struct associated to the path

function [ dPath, phiPath, splineData ] = loadPathDiff(filename, varargin)

autocomplete = true;
workdir = '';

% Some code for handling optional inputs
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'autocomplete'
                ii = ii + 1;
                if isa(varargin{ii}, 'integer') || ...
                        isa(varargin{ii}, 'logical')
                    autocomplete = varargin{ii};
                else
                    error('Invalid value for option ''autocomplete''.');
                end
            case 'workdir'
                ii = ii + 1;
                if isa(varargin{ii}, 'char')
                    workdir = varargin{ii};
                else
                    error('Invalid value for option ''workdir''.');
                end
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end
    ii = ii + 1; 
    end
end

if ~isempty(workdir) && workdir(end) ~= '/'
    workdir(end+1) = '/';
end
filename = [workdir, filename];

if exist(filename, 'file') % Try direct existence first
    curveFile = matfile(filename);
elseif autocomplete % Look for files matching pattern
    candidates = dir([filename, '*.mat']);
    if ~isempty(candidates)
        curveFile = matfile([workdir, candidates(1).name]);
    else
        error('File ''%s'' not found.', filename);
    end
else
    error('File ''%s'' not found.', filename);
end

if ~isempty(who(curveFile, 'dPath'))
    dPath = curveFile.dPath;
else
    error('File ''%s'' does not contain field d0.', filename);
end

if ~isempty(who(curveFile, 'phiPath'))
    phiPath = curveFile.phiPath;
else
    error('File ''%s'' does not contain field d0.', filename);
end

if nargout > 2
    if ~isempty(who(curveFile, 'splineData'))
        splineData = curveFile.splineData;
    else
        error('File ''%s'' does not contain field splineData.', filename);
    end
end

end
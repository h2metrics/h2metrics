%% savePathDiff
% Saves a path, a path of diffeomorphisms and associated splineData to a 
% file.
%
% The filename is constructed as
%       name_n4_N20_nT2_Nt10_diff.mat
%
% Input
%   name
%       Path name used to construct filename by appending info about
%       N, nS, Nt and nT.
%   dPath
%       Path to save
%   phiPath
%       Diffeomorphisms to save
%   splineData
%       splineData associated to the path
%
% Optional parameters
%   'workdir'
%       Appends an optional working directory in front of filename

function savePathDiff( name, dPath, phiPath, splineData, varargin )

workdir = '';

% Some code for handling optional inputs
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
         switch (lower(varargin{ii}))
             case 'workdir'
                 ii = ii + 1;
                 if isa(varargin{ii}, 'char')
                     workdir = varargin{ii};
                 else
                     error('Invalid value for option ''workDir''');
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

N = splineData.N;
Nt = splineData.Nt;
nS = splineData.nS;
nT = splineData.nT;

filename = [workdir, name,'_n',num2str(nS),'_N',num2str(N),...
            '_nT',num2str(nT),'_Nt',num2str(Nt),'_diff.mat'];

save(filename, 'dPath', 'phiPath', 'splineData');
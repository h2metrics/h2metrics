%% saveCurve
% Saves a curve and associated splineData to a file.
%
% The filename is constructed as
%       name_n4_N20.mat
%
%
% Input
%   name
%       Curve name used to construct filename by appending info about
%       N and nS
%   d0
%       Curve to save
%   splineData
%       splineData associated to curve
%
% Optional parameters
%   'workdir'
%       Appends an optional working directory in front of filename
%   'manualName' = {true, false}
%       Don't append anything to given filename

function saveCurve( name, d0, splineData, varargin )

workdir = '';
manualName = false;

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
                     error('Invalid value for option ''workDir''.');
                 end
             case 'manualname'
                 ii = ii + 1;
                 if isa(varargin{ii}, 'integer') || ...
                     isa(varargin{ii}, 'logical')
                     manualName = logical(varargin{ii});
                 else
                     error('Invalid value for option ''manualName''.');
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
nS = splineData.nS;

if manualName
    filename = [workdir, name];
else
    filename = [workdir, name,'_n',num2str(nS),'_N',num2str(N),'.mat'];
end

save(filename, 'd0', 'splineData');

end
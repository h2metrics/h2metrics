%% loadDataSet
% Function loads the spline approximation for a given dataset
%
% Input
%   setName
%       Describes the dataset to be loaded. At the moment the following are
%       supported:
%           hela_murphy, xavier_heart, corpus_callosum_tilak
%   splineData
%       Describes the spline degree to be used.
%   dataDir
%       Location of the source data.
%
% Optional arguments
%       The precise set of optional arguments is different for each
%       dataset. Please refer to the documentation of the function
%       loadDataSetXYZ for details. Some common parameters will be
%       explained here.
%   'reloadData'
%       Forces the recreation of splines from the source data. Existing
%       splines are not used.
%   'noPlot'
%       Skips the plotting of the data. Can be used to make code faster.
%   'constSpeed'
%       Parametrizes spline to constant speed.
%   'ind' = Array
%       Only returns the data corresponding to the given indices. The
%       ordering of the curves will depend on the dataset. Usually it will
%       correspond to the alphabetic ordering of files in the source
%       directory. Example
%           loadDataSet( ..., 'ind', [1, 2, 5, 6])
%       To return all curves, pass an empty array ([]).
%
% Output
%   dList
%       Cell array with the spline curves.
function dList = loadDataSet( setName, splineData, dataDir, varargin )

switch (lower(setName))
    case 'hela_murphy'
        dList = loadDataSetHela(splineData, dataDir, varargin{:});
    case 'xavier_heart'
        dList = loadDataSetXavier(splineData, dataDir, varargin{:});
    case 'corpus_callosum_tilak'
        dList = loadDataSetTilak(splineData, dataDir, varargin{:});
    case 'basic'
        dList = loadDataSetBasic(splineData, varargin{:});
    otherwise
        error('Invalid dataset: ''%s''.',varargin{ii});
end
    
end





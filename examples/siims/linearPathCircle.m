%% linearPath
%
% Function computes the linear interpolation between d0 and d1 using the
% spline described by knotsT, Nt and nT.
%
% Input
%   d0, d1
%       Matrices to be linearly interpolated
%   splineData
%       Contains knotsT, Nt and nT
%
% Output
%   dPath
%       The interpolated path. The dimensions of dPath are
%           [size(d0,1)*Nt, size(d0,2)]

function dPath = linearPathCircle( d0, d1,splineData )

if ~isequal( size(d0), size(d1) )
    error('Dimension mismatch.');
end
options = struct( 'optDiff', true, ...
                  'optTra', true, ...
                  'optRot', true, ...
                  'optShift', true, ...
                  'tolFun', 1e-12, ...
                  'tolX', 1e-12, ...
                  'display', 'iter-detailed', ... % 'iter-detailed'
                  'maxIter', 400 );

options.rigidA = [1, 0, 0];    
              
              
N = size(d0, 1);
radius1 = max(max(d0(:,1))-min(d0(:,1)),max(d0(:,2))-min(d0(:,2)));
radius2 =max(max(d1(:,1))-min(d1(:,1)),max(d1(:,2))-min(d1(:,2)));
radius=max(radius1,radius2)/2;


fmiddle = @(t) radius*[sin(t), cos(t)] ;
dmiddle = constructSplineApproximation(fmiddle,splineData);
[dAligned, ~] = rigidAlignment( {d0,dmiddle}, splineData,'options', options); 
[dAligned2, ~] = rigidAlignment( {dAligned{2},d1}, splineData,'options', options);
                                         
splineData1 = constructEmptySplineData;
splineData1.N = splineData1.N ;%no. control points, must be bigger than n+1
splineData1.Nt = splineData.Nt/2; %Number of time control points
splineData1.Nphi = splineData.Nphi; %No. control points for diffeomorphisms
splineData1.nS = splineData.nS; %spacial degree
splineData1.nT = splineData.nT; %time degree
splineData1.nPhi = splineData.nPhi; %diffemorphism degree
splineData1.quadDegree = splineData.quadDegree; %Quadrature precission
splineData1.dSpace = splineData.dSpace;
splineData1.noInterpolS = splineData.noInterpolS ; % For composition
splineData1 = constructKnots(splineData1);  

dPath1 = linearPath( dAligned{1}, dAligned{2}, splineData1);
dPath2 = linearPath( dAligned2{1}, dAligned2{2}, splineData1);
dPath=[dPath1; dPath2];
end




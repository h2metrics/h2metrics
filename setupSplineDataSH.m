function splineData = setupSplineDataSH(task)
%setupSplineDataSH(task) sets up the collocation and weight matrices used for plotting. 
%
% TODO: Summary

splineData = struct(...
    'KT',   [],...  % time knots
    'KU',   [],...  % spacial knots
    'innerKT', [],...
    'innerKU', [],...
    'QT',   [],...  % time quadrature points
    'QU',   [],...  % spacial quadrature points
    'nQT',  [],...  % number of time quadrature points
    'nQU',  [],...  % number of spacial quadrature points
    'WT',   [],...  % time quadrature weights
    'WU',   [],...  % spacial quadrature weights
    'W',    [],...  % outer product of the two weight vectors
    'B',    [],...  % evaluation of tensor product splines at quadrature points
    'Bt',   [],...  % one time derivative
    'Btu',  [],...  % one time, one space derivative
    'Btuu', [],...  % etc
    'Bu',   [],...
    'Buu',  [],...
    'Buuu', [],...
    'BU',   [],...  % only the spatial splines
    'BUu',  [],...  % one derivative
    'BUuu', [],...  % two derivatives
    'BUuuu',[]...  % etc
);

% calculate time and space knots
splineData.KT = task.splineData.KT;
splineData.KU = task.splineData.KU;
splineData.innerKT = task.splineData.innerKT;
splineData.innerKU = task.splineData.innerKU;

% calculate quadrature points and weights
splineData.QU = linspace(0,2*pi,task.nrPlotPoints(1));
splineData.QT = linspace(0,1,task.nrPlotPoints(2));
splineData.nQU = length(splineData.QU);
splineData.nQT = length(splineData.QT);

splineData.WU = ones([1,task.nrPlotPoints(1)]);
splineData.WU(1) = 0.5; splineData.WU(end) = 0.5;
splineData.WT = ones([1,task.nrPlotPoints(2)]);
splineData.WT(1) = 0.5; splineData.WT(end) = 0.5;
splineData.W = (splineData.WU'*splineData.WT)'; 

% calculate 0,1,2,3 derivatives of spline basis in space
B_S_quad = spcol( splineData.KU, task.dU+1, brk2knt( splineData.QU, 4 ),'sparse');

%Periodic B-splines
B_S_quad_per = [B_S_quad(:,1:task.dU) + B_S_quad(:,end-task.dU+1:end), B_S_quad(:,task.dU+1:end-task.dU)];

% calculate 0,1 derivatives of spline basis in time
B_T_quad = spcol( splineData.KT, task.dT+1, brk2knt( splineData.QT, 2 ),'sparse');

% extract the spatial collocation matrices
splineData.BU    = B_S_quad_per(1:4:end,:);
splineData.BUu   = B_S_quad_per(2:4:end,:);
splineData.BUuu  = B_S_quad_per(3:4:end,:);
splineData.BUuuu = B_S_quad_per(4:4:end,:);

% calculate the tensor product spline matrices
for ii = task.nKU:-1:1
    for jj = task.nKT:-1:1;
     
    splineData.B(:,ii+(jj-1)*task.nKU) = reshape(B_S_quad_per(1:4:end,ii)*B_T_quad(1:2:end,jj)',...
        splineData.nQU*splineData.nQT ,1) ;
    splineData.Bt(:,ii+(jj-1)*task.nKU) = reshape(B_S_quad_per(1:4:end,ii)*B_T_quad(2:2:end,jj)',...
        splineData.nQU*splineData.nQT ,1)' ;
    splineData.Btu(:,ii+(jj-1)*task.nKU) = reshape(B_S_quad_per(2:4:end,ii)*B_T_quad(2:2:end,jj)',...
        splineData.nQU*splineData.nQT ,1)' ;
    splineData.Btuu(:,ii+(jj-1)*task.nKU) = reshape(B_S_quad_per(3:4:end,ii)*B_T_quad(2:2:end,jj)',...
        splineData.nQU*splineData.nQT ,1)' ;
    splineData.Bu(:,ii+(jj-1)*task.nKU) = reshape(B_S_quad_per(2:4:end,ii)*B_T_quad(1:2:end,jj)',...
        splineData.nQU*splineData.nQT ,1)' ;
    splineData.Buu(:,ii+(jj-1)*task.nKU) = reshape(B_S_quad_per(3:4:end,ii)*B_T_quad(1:2:end,jj)',...
        splineData.nQU*splineData.nQT ,1)' ;
    splineData.Buuu(:,ii+(jj-1)*task.nKU) = reshape(B_S_quad_per(4:4:end,ii)*B_T_quad(1:2:end,jj)',...
        splineData.nQU*splineData.nQT ,1)' ;
    end
end

end




clear all
close all

% Load both curves. Each is a structure containing a field 'x' with the
% vertex positions and a field 'G' giving the connetivity matrix.
load curve3.mat 
load curve4.mat

% Parameters for the kernels in the metric
objfun.kernel_geom='gaussian';
objfun.kernel_size_geom=20;
objfun.kernel_grass='binet';

% Varifold distance
tic
g=varifoldnorm(fs3,fs4,objfun)
toc

% Gradient of the distance wrt vertices positions
tic
[dxg]= dvarifoldnorm(fs3,fs4,objfun);
toc

figure()
hold on
plot(fs3.x(:,1),fs3.x(:,2),'b')
plot(fs4.x(:,1),fs4.x(:,2),'r')
quiver(fs3.x(:,1),fs3.x(:,2),-dxg(:,1),-dxg(:,2),'k')
axis equal
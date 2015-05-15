%% Mutlidimensional scaling - eigenvectors, dist
disp('Riemannian distance MDS, eigenvectors');

load('DistA1B01C0.mat');
tmpDist = Dist;

[Y, e] = cmdscale(tmpDist);

%% Multidimensional scaling - 2d plot, length
close all;
axis off;
hold on;
disp('Riemannian distance MDS, 2d plot');
n1=9;
n2=length(Y(:,1));
deltax=0.02;
deltay=-0.007;
plot(Y(1:n1,1), Y(1:n1,2), 'o','Color','blue');
plot(Y(n1+1:n2,1), Y(n1+1:n2,2), 'x','Color','red');
hold on;
for kk = 1:n1
    text(Y(kk,1)+deltax, Y(kk,2)+deltay, num2str(kk),'Color','blue');
end
for kk = n1+1:n2
    text(Y(kk,1)+deltax, Y(kk,2)+deltay, num2str(kk),'Color','red');
end
set(gcf,'color','w');
hold off;




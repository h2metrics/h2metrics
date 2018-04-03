%% Cluster and dentogramm
load('DistA1B01Ce-9.mat');
Y = squareform(Dist,'tovector');
Z = linkage(Y);
c = cluster(Z,'maxclust',2)
dendrogram(Z)

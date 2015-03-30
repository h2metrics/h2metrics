function writedatfile(task)
%writedatfile writes the .dat file needed by AMPL
%   writedatfile(d0, d1, quadData, datfile) creates a .dat file 
%   with the name filename. d0 and d1 are the initial and final controls, 
%   respectively, and quadData is the quadrature data structure.

%% Determine the number of control and quadrature points
nKT = task.nKT;
nKU = task.nKU;
nQT = task.splineData.nQT;
nQU = task.splineData.nQU;

%% Open the .dat file
fileid=fopen(task.datfile,'w');

%% Write the scalar parameters
fprintf(fileid,'param pHor := %+24.17e; \n',task.pHor);
fprintf(fileid,'param pH0 := %+24.17e; \n',task.pH0);
fprintf(fileid,'param pKappa := %+24.17e; \n',task.pKappa);
fprintf(fileid,'param pH2 := %+24.17e; \n',task.pH2);
fprintf(fileid,'param nKT := %d; # number of t-knots\n',nKT);
fprintf(fileid,'param nKU := %d; # number of p-knots\n',nKU);
fprintf(fileid,'param nQT := %d; # number of t-quadrature points\n',nQT);
fprintf(fileid,'param nQU := %d; # number of p-quadrature points\n',nQU);
fprintf(fileid,'\n');

%% Write the initial and final controls d0, d1
fprintf(fileid,'param d0: 1 2 := \n');
    for i=1:nKU 
        fprintf(fileid,'%d\t%+24.17e\t%+24.17e\n',i,task.results.d0(i,1),task.results.d0(i,2));
    end
fprintf(fileid,';\n\n');

fprintf(fileid,'param d1: 1 2 := \n');
    for i=1:nKU 
        fprintf(fileid,'%d\t%+24.17e\t%+24.17e\n',i,task.results.d1(i,1),task.results.d1(i,2));
    end
fprintf(fileid,';\n\n');

%% Write the spline collocation matrices
fprintf(fileid,'\nparam: iB: B :=\n');
indices = find(task.splineData.B);
values = nonzeros(task.splineData.B);
[u,t,uu,tt] = ind2sub([nQU, nQT, nKU, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%+24.17e\n',[t,u,tt,uu,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBt: Bt :=\n');
indices = find(task.splineData.Bt);
values = nonzeros(task.splineData.Bt);
[u,t,uu,tt] = ind2sub([nQU, nQT, nKU, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%+24.17e\n',[t,u,tt,uu,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBtu: Btu :=\n');
indices = find(task.splineData.Btu);
values = nonzeros(task.splineData.Btu);
[u,t,uu,tt] = ind2sub([nQU, nQT, nKU, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%+24.17e\n',[t,u,tt,uu,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBtuu: Btuu :=\n');
indices = find(task.splineData.Btuu);
values = nonzeros(task.splineData.Btuu);
[u,t,uu,tt] = ind2sub([nQU, nQT, nKU, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%+24.17e\n',[t,u,tt,uu,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBu: Bu :=\n');
indices = find(task.splineData.Bu);
values = nonzeros(task.splineData.Bu);
[u,t,uu,tt] = ind2sub([nQU, nQT, nKU, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%+24.17e\n',[t,u,tt,uu,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBuu: Buu :=\n');
indices = find(task.splineData.Buu);
values = nonzeros(task.splineData.Buu);
[u,t,uu,tt] = ind2sub([nQU, nQT, nKU, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%+24.17e\n',[t,u,tt,uu,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBuuu: Buuu :=\n');
indices = find(task.splineData.Buuu);
values = nonzeros(task.splineData.Buuu);
[u,t,uu,tt] = ind2sub([nQU, nQT, nKU, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%+24.17e\n',[t,u,tt,uu,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBU: BU :=\n');
indices = find(task.splineData.BU);
values = nonzeros(task.splineData.BU);
[u,uu] = ind2sub([nQU, nKU], indices);
fprintf(fileid,'%d\t%d\t%+24.17e\n',[u,uu,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBUu: BUu :=\n');
indices = find(task.splineData.BUu);
values = nonzeros(task.splineData.BUu);
[u,uu] = ind2sub([nQU, nKU], indices);
fprintf(fileid,'%d\t%d\t%+24.17e\n',[u,uu,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBUuu: BUuu :=\n');
indices = find(task.splineData.BUuu);
values = nonzeros(task.splineData.BUuu);
[u,uu] = ind2sub([nQU, nKU], indices);
fprintf(fileid,'%d\t%d\t%+24.17e\n',[u,uu,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBUuuu: BUuuu :=\n');
indices = find(task.splineData.BUuuu);
values = nonzeros(task.splineData.BUuuu);
[u,uu] = ind2sub([nQU, nKU], indices);
fprintf(fileid,'%d\t%d\t%+24.17e\n',[u,uu,values]');
fprintf(fileid,';\n\n');


%% Write the weight matrix

fprintf(fileid,'\nparam W :=\n');
indices = find(task.splineData.W);
values = nonzeros(task.splineData.W);
[u,t] = ind2sub([nQU, nQT], indices);
fprintf(fileid,'%d\t%d\t%+24.17e\n',[t,u,values]');
fprintf(fileid,';\n\n');

%% close the file
fclose(fileid);
end

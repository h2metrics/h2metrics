function writedatfile(d0,d1,quadData,datfile)
%writedatfile writes the .dat file needed by AMPL
%   writedatfile(d0, d1, quadData, datfile) creates a .dat file 
%   with the name filename. d0 and d1 are the initial and final controls, 
%   respectively, and quadData is the quadrature data structure.

%% Determine the number of control and quadrature points
nKP = quadData.noControlPoints(1);
nKT = quadData.noControlPoints(2);
nQP = length(quadData.points{1});
nQT = length(quadData.points{2});

%% Open the .dat file
fileid=fopen(datfile,'w');

%% Write the scalar parameters
fprintf(fileid,'param nKP := %d; # number of p-knots\n',nKP);
fprintf(fileid,'param nKT := %d; # number of t-knots\n',nKT);
fprintf(fileid,'param nQP := %d; # number of p-quadrature points\n',nQP);
fprintf(fileid,'param nQT := %d; # number of t-quadrature points\n',nQT);
fprintf(fileid,'\n');

%% Write the initial and final controls d0, d1
fprintf(fileid,'param d0: 1 2 := \n');
    for i=1:nKP 
        fprintf(fileid,'%d\t%.17e\t%.17e\n',i,d0(i,1),d0(i,2));
    end
fprintf(fileid,';\n\n');

fprintf(fileid,'param d1: 1 2 := \n');
    for i=1:nKP 
        fprintf(fileid,'%d\t%.17e\t%.17e\n',i,d1(i,1),d1(i,2));
    end
fprintf(fileid,';\n\n');

%% Write the spline collocation matrices
fprintf(fileid,'\nparam: iB: B :=\n');
indices = find(quadData.B);
values = nonzeros(quadData.B);
[p,t,pp,tt] = ind2sub([nQP, nQT, nKP, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%.17e\n',[t,p,tt,pp,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBp: Bp :=\n');
indices = find(quadData.Bu);
values = nonzeros(quadData.Bu);
[p,t,pp,tt] = ind2sub([nQP, nQT, nKP, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%.17e\n',[t,p,tt,pp,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBt: Bt :=\n');
indices = find(quadData.Bt);
values = nonzeros(quadData.Bt);
[p,t,pp,tt] = ind2sub([nQP, nQT, nKP, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%.17e\n',[t,p,tt,pp,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBpt: Bpt :=\n');
indices = find(quadData.But);
values = nonzeros(quadData.But);
[p,t,pp,tt] = ind2sub([nQP, nQT, nKP, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%.17e\n',[t,p,tt,pp,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBpp: Bpp :=\n');
indices = find(quadData.Buu);
values = nonzeros(quadData.Buu);
[p,t,pp,tt] = ind2sub([nQP, nQT, nKP, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%.17e\n',[t,p,tt,pp,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBppt: Bppt :=\n');
indices = find(quadData.Buut);
values = nonzeros(quadData.Buut);
[p,t,pp,tt] = ind2sub([nQP, nQT, nKP, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%.17e\n',[t,p,tt,pp,values]');
fprintf(fileid,';\n\n');


%% Write the weight matrix

fprintf(fileid,'\nparam W :=\n');
quadWeightMult = (quadData.weights{2}'*quadData.weights{1})'; %Multiply u t(i)
indices = find(quadWeightMult);
values = nonzeros(quadWeightMult);
[p,t] = ind2sub([nQP, nQT], indices);
fprintf(fileid,'%d\t%d\t%.17e\n',[t,p,values]');
fprintf(fileid,';\n\n');

%% close the file
fclose(fileid);
end

function writedatfile(d0,d1,splineData,quadData,quadDataTensor,datfile)
%writedatfile writes the .dat file needed by AMPL
%   writedatfile(d0, d1, quadData, datfile) creates a .dat file 
%   with the name filename. d0 and d1 are the initial and final controls, 
%   respectively, and quadData and quadDataTensor are the quadrature data structure.

%% Determine the number of control and quadrature points
nKP = splineData.N;
nKT = splineData.Nt;
nQP =  length(quadData.quadPointsS);   
nQT =  length(quadData.quadPointsT);

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
indices = find(quadDataTensor.B);
values = nonzeros(quadDataTensor.B);
[p,t,pp,tt] = ind2sub([nQP, nQT, nKP, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%.17e\n',[t,p,tt,pp,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBp: Bp :=\n');
indices = find(quadDataTensor.Bu);
values = nonzeros(quadDataTensor.Bu);
[p,t,pp,tt] = ind2sub([nQP, nQT, nKP, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%.17e\n',[t,p,tt,pp,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBt: Bt :=\n');
indices = find(quadDataTensor.Bt);
values = nonzeros(quadDataTensor.Bt);
[p,t,pp,tt] = ind2sub([nQP, nQT, nKP, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%.17e\n',[t,p,tt,pp,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBpt: Bpt :=\n');
indices = find(quadDataTensor.But);
values = nonzeros(quadDataTensor.But);
[p,t,pp,tt] = ind2sub([nQP, nQT, nKP, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%.17e\n',[t,p,tt,pp,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBpp: Bpp :=\n');
indices = find(quadDataTensor.Buu);
values = nonzeros(quadDataTensor.Buu);
[p,t,pp,tt] = ind2sub([nQP, nQT, nKP, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%.17e\n',[t,p,tt,pp,values]');
fprintf(fileid,';\n\n');

fprintf(fileid,'\nparam: iBppt: Bppt :=\n');
indices = find(quadDataTensor.Buut);
values = nonzeros(quadDataTensor.Buut);
[p,t,pp,tt] = ind2sub([nQP, nQT, nKP, nKT], indices);
fprintf(fileid,'%d\t%d\t%d\t%d\t%.17e\n',[t,p,tt,pp,values]');
fprintf(fileid,';\n\n');


%% Write the weight matrix

fprintf(fileid,'\nparam W :=\n');
quadWeightMult = (quadData.quadWeightsS*quadData.quadWeightsT'); %Multiply u t(i)
indices = find(quadWeightMult);
values = nonzeros(quadWeightMult);
[p,t] = ind2sub([nQP, nQT], indices);
fprintf(fileid,'%d\t%d\t%.17e\n',[t,p,values]');
fprintf(fileid,';\n\n');

%% close the file
fclose(fileid);
end

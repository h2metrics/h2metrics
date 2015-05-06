function writedatfile2(splineData,quadData,quadDataTensor,datfile)
%writedatfile writes the .dat file containing the collocation matrices.
%   writedatfile(splineData,quadData,quadDataTensor, datfile) creates a .dat file 
%   with the name filename.
%   splineData, quadData and quadDataTensor are the quadrature data structure.

%% Determine the number of control and quadrature points
nKP = splineData.N;
nKT = splineData.Nt;
nQP =  length(quadData.quadPointsS);   
nQT =  length(quadData.quadPointsT);


%% Open the .dat file
fileid=fopen(datfile,'w');


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


%% close the file
fclose(fileid);
end

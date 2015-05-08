function writerunfile(runfile, modfile, datfile1,datfile2, tabfile, tabledef, options)
%writerunfile writes the .run file needed by AMPL
%   writedatfile(runfile, modfile, datfile, tabfile, tabledef, options) 
%   runfile etc are the names of the .run, .mod, .dat, and .tab files.
%   tabledef is a cell array containing the AMPL code and column headers of the table definition
%   options is a cell array containing options passed to AMPL and/or Ipopt
fileid=fopen(runfile,'w');
fprintf(fileid,'model %s;\n', modfile);
fprintf(fileid,'data %s;\n', datfile1);
fprintf(fileid,'data %s;\n', datfile2);
fprintf(fileid,'option solver ipopt;\n');
for i=1:length(options)
  fprintf(fileid,'%s;\n',options{i});
end
[tabfilepath,tabfilename,tabfileext]=fileparts(tabfile);
fprintf(fileid,'table %s OUT:{t in KT,p in KP} -> \n  [KT,KP],\n',tabfilename);
fprintf(fileid,'  %s ~ %s',tabledef{1,1},tabledef{1,2});
for i=2:length(tabledef)
  fprintf(fileid,',\n');
  fprintf(fileid,'  %s ~ %s',tabledef{i,1},tabledef{i,2});
end
fprintf(fileid,'\n;\nsolve;\nwrite table %s;\n\n',tabfilename); 
fprintf(fileid,'display(energy);\n\n');

fclose(fileid);
end 

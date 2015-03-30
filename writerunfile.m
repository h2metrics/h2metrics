function writerunfile(task)
%writerunfile writes the .run file needed by AMPL
%   writedatfile(runfile, modfile, datfile, tabfile, tabledef, options) 
%   runfile etc are the names of the .run, .mod, .dat, and .tab files.
%   tabledef is a cell array containing the AMPL code and column headers of the table definition
%   options is a cell array containing options passed to AMPL and/or Ipopt
fileid=fopen(task.runfile,'w');
fprintf(fileid,'model %s;\n', task.modfile);
fprintf(fileid,'data %s;\n', task.datfile);
fprintf(fileid,'option solver ipopt;\n');
for i=1:length(task.amploptions)
  fprintf(fileid,'%s;\n',task.amploptions{i});
end
[tabfilepath,tabfilename,tabfileext]=fileparts(task.tabfile);
fprintf(fileid,'table %s OUT:{t in KT,u in KU} -> \n  [KT,KU],\n',tabfilename);
fprintf(fileid,'  %s ~ %s',task.tabledef{1,1},task.tabledef{1,2});
for i=2:length(task.tabledef)
  fprintf(fileid,',\n');
  fprintf(fileid,'  %s ~ %s',task.tabledef{i,1},task.tabledef{i,2});
end
fprintf(fileid,'\n;\nsolve;\nwrite table %s;\n\n',tabfilename); 
fprintf(fileid,'display(energy);\n\n');

fclose(fileid);
end 

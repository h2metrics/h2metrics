function [data, variables] = readtabfile(filename)
%readtabfile Read .tab file from AMPL
%   [data, variables] = readtabfile(filename) reads the specified .tab file
%   into the Matlab array 'data' with size(data) = [len(variables), N_V, N_S+1]. 
%   The column headers of the .tab file are stored in the vector 'variables'. 

fileid=fopen(filename);
% set n(1) to the number of integer and n(2) to the number of double columns
n=fscanf(fileid, 'ampl.tab %d %d ');
% read the variable names from line 2 of the .tab file
l1=fgetl(fileid);
variables=strread(l1,'%s');
% discard the parameter columns
variables = variables(n(1)+1:end);

% create format string to read subsequent lines
format='';
for i = 1:n(1)
    format = [format, '%d '];
end
for i = 1:n(2)
    format = [format, '%f '];
end

% read data, close file
data=fscanf(fileid,format,[n(1)+n(2),Inf]);
fclose(fileid);

% get the number of vertices N(1), and time points N(2) 
for i = 1: n(1)
    N(i) = max(data(n(1)+1-i,:));
end

% discard the parameter columns
data = data(n(1)+1:end,:);

% check if c has the correct size
if size(data) ~= [n(2), prod(N)]
    error('Too many or not enough lines found in the .tab file');
end
data=reshape(data,[n(2) N]);
end

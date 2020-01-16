function bvals=read_bval(file_loc)
%%
% reads b-values from text file
%
% Code is written by Oliver Gurney-Champion
% o.j.gurney-chapion@amsterdamumc.nl
%%
txts=fopen(file_loc);
tline = fgetl(txts);
bvals=str2num(tline);
fclose(txts);

end
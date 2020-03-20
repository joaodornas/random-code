
filename = 'myfile.txt';
fid = fopen(filename, 'wt');
disp(num2str(fid));
nbytes = fprintf(fid,'%s','isso eh um teste');
fclose(fid);
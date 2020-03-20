
N = nifti;
N.dat = file_array(fname,dim,dtype,offset,scl_slope,scl_inter);
N.mat = nifti_file.mat;
N.mat_intent = nifti_file.mat_intent;
N.mat0 = nifti_file.mat0;
N.mat0_intent = nifti_file.mat0_intent;
N.descrip = descrip;
create(N);
dat = N.dat;
dat.scl_slope = scl_slope;
dat.scl_inter = scl_inter;

if length(size(input_data)) == 2
    
    dat(:,:) = input_data;

elseif length(size(input_data)) == 3
    
    dat(:,:,:) = input_data;
    
elseif length(size(input_data)) == 4
    
    dat(:,:,:,:) = input_data;
    N.timing = nifti_file.timing;
    
end

clear N
clear dat
clear input_data
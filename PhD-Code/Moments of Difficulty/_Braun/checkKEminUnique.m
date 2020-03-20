
start_ = input('Start:');

end_ = input('End:');

prefix = 'MOT-Analyze-V2-10-m-';

for i=start_:end_
    i
    datafile = strcat(prefix,int2str(i),'.mat');
    
    info = load(datafile);
    
    window.kentropyColorMin(i) = info.kECmin;
    window.kentropyMotionMin(i) = info.kEMmin;
    window.kentropyBounceMin(i) = info.kEBmin;
    window.kentropyShiftMin(i) = info.kEQmin;
    
    window.kweightColorMin(i) = info.kVCmin;
    window.kweightMotionMin(i) = info.kVMmin;
    window.kweightBounceMin(i) = info.kVBmin;
    window.kweightShiftMin(i) = info.kVQmin;
    
    window.kentropyColorMax(i) = info.kECmax;
    window.kentropyMotionMax(i) = info.kEMmax;
    window.kentropyBounceMax(i) = info.kEBmax;
    window.kentropyShiftMax(i) = info.kEQmax;
    
    window.kweightColorMax(i) = info.kVCmax;
    window.kweightMotionMax(i) = info.kVMmax;
    window.kweightBounceMax(i) = info.kVBmax;
    window.kweightShiftMax(i) = info.kVQmax;
    
    clear info
    clear datafile
   
end

[uniECmin, nuniECmin] = count_unique(window.kentropyColorMin);
[uniEMmin, nuniEMmin] = count_unique(window.kentropyMotionMin);

[uniVCmin, nuniVCmin] = count_unique(window.kweightColorMin);
[uniVMmin, nuniVMmin] = count_unique(window.kweightMotionMin);

[uniECmax, nuniECmax] = count_unique(window.kentropyColorMax);
[uniEMmax, nuniEMmax] = count_unique(window.kentropyMotionMax);

[uniVCmax, nuniVCmax] = count_unique(window.kweightColorMax);
[uniVMmax, nuniVMmax] = count_unique(window.kweightMotionMax);


idx_unique_ECmin = find(nuniECmin == 1);
k_unique_ECmin = uniECmin(idx_unique_ECmin);
idx_MOT_unique_ECmin = find(ismember(window.kentropyColorMin,k_unique_ECmin));

idx_unique_EMmin = find(nuniEMmin == 1);
k_unique_EMmin = uniEMmin(idx_unique_EMmin);
idx_MOT_unique_EMmin = find(ismember(window.kentropyMotionMin,k_unique_EMmin));

idx_unique_ECmax = find(nuniECmax == 1);
k_unique_ECmax = uniECmax(idx_unique_ECmax);
idx_MOT_unique_ECmax = find(ismember(window.kentropyColorMax,k_unique_ECmax));

idx_unique_EMmax = find(nuniEMmax == 1);
k_unique_EMmax = uniEMmax(idx_unique_EMmax);
idx_MOT_unique_EMmax = find(ismember(window.kentropyMotionMax,k_unique_EMmax));


idx_unique_VCmin = find(nuniVCmin == 1);
k_unique_VCmin = uniVCmin(idx_unique_VCmin);
idx_MOT_unique_VCmin = find(ismember(window.kweightColorMin,k_unique_VCmin));

idx_unique_VMmin = find(nuniVMmin == 1);
k_unique_VMmin = uniVMmin(idx_unique_VMmin);
idx_MOT_unique_VMmin = find(ismember(window.kweightMotionMin,k_unique_VMmin));

idx_unique_VCmax = find(nuniVCmax == 1);
k_unique_VCmax = uniVCmax(idx_unique_VCmax);
idx_MOT_unique_VCmax = find(ismember(window.kweightColorMax,k_unique_VCmax));

idx_unique_VMmax = find(nuniVMmax == 1);
k_unique_VMmax = uniVMmax(idx_unique_VMmax);
idx_MOT_unique_VMmax = find(ismember(window.kweightMotionMax,k_unique_VMmax));

save('uniqueKEmin.mat');
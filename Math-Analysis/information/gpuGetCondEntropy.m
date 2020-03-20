function [ HCond, HCondMAbound ] = gpuGetCondEntropy(CondFrequencies)

for i=1:length(CondFrequencies)

    frequencies = CondFrequencies(i).freq;

    [ H, HMAbound ] = gpuGetEntropy(frequencies);
    
    HCond(i) = H;
    
    HCondMAbound(i) = HMAbound;
    
    clear H;
    clear HMAbound;
    
end


end


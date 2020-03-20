function [ H, HMAbound ] = gpuGetEntropy(frequencies)

H = 0;

for f=1:length(frequencies)

   H = H - frequencies(f)*log2(frequencies(f));

end

sumOfFrequencies = 0;
for f=1:length(frequencies)

   sumOfFrequencies = sumOfFrequencies + frequencies(f)^2;

end

HMAbound = - log2(sumOfFrequencies);


end


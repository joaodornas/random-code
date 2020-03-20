function TB_optimal = TrialsFFoptimal(TB)

time_bins = size(TB,1);
nTrials = size(TB,2);

for i=1:time_bins

    for k=1:nTrials

       for j=1:i

           nBins = j;

           TB_optimal(i).Trial(k).bin(j) = TB(i,k,j);

       end

    end

end


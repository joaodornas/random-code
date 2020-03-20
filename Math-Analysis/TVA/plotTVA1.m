function plotTVA1

ED1 = [10	20	30	40	50	100	120	200];

Score0 = [99	45	14	6	7	3	4	6];

Score1 = [1	55	86	94	93	97	96	94];

plot(ED1(1:end),Score0(1:end),'bo-');

hold on

plot(ED1(1:end),Score1(1:end),'ro-');

xlabel('Exposure Duration (ms)');

ylabel('Probability of Reporting');


end


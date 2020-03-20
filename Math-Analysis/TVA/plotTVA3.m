function plotTVA3

ED1 = [10	20	30	40	50	100	120	200];

Score0 = [90	27	3	4	1	1	0	0];

Score1 = [10	58	42	25	16	5	7	1];

Score2 = [0 13 48 51 42 23 21 24];

Score3 = [0 2 7 20 41 71 72 75];

plot(ED1(1:end),Score0(1:end),'b-o');

hold on

plot(ED1(1:end),Score1(1:end),'r-o');

hold on

plot(ED1(1:end),Score2(1:end),'g-o');

hold on

plot(ED1(1:end),Score3(1:end),'k-o');


end


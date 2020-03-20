function plotTVA2

ED1 = [10	20	30	40	50	100	120	200];

Score0 = [93	28	4	4	1	0	2	0];

Score1 = [7	62	44	24	17	11	7	9];

Score2 = [0 10 52 72 82 89 91 91];

plot(ED1(1:end),Score0(1:end),'b-o');

hold on

plot(ED1(1:end),Score1(1:end),'r-o');

hold on

plot(ED1(1:end),Score2(1:end),'g-o');


end


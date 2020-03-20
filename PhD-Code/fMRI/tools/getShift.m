function [shift_pos, shift_neg] = getShift(dist1,dist2)

step = 0.001;
interval = -1:step:1;

mid = round(length(interval)/2);

pos1 = dist1(mid:end);
int_pos1 = sum(pos1);
pos2 = dist2(mid:end);
int_pos2 = sum(pos2);

neg1 = dist1(1:mid);
int_neg1 = sum(neg1);
neg2 = dist2(1:mid);
int_neg2 = sum(neg2);

shift_pos = int_pos1-int_pos2;
shift_neg = int_neg1-int_neg2;

end
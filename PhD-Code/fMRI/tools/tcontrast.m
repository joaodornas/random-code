function [H, P, threshold_labels] = tcontrast(ROIs1, ROIs2, labels, alpha)

contrast = ROIs1 - ROIs2;

[H, P] = ttest(contrast);

idx = find(P < alpha);

threshold_labels = labels(idx);

end


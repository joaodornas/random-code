

katha = [3 1 2 2 ]
robert = [2 3 2 2]
jan = [3 1 2 3]
julia = [1 3 1]

all_subjects = [katha,robert,jan,julia]

mean(all_subjects)

mn = 2

all_subjects_zs = all_subjects - mn

[H,P] = ttest(all_subjects_zs)




%%% Columns

%% column 1 - Passive - All - Down
%% column 2 - Passive - All - Up
%% column 3 - Passive - Neg - Down
%% column 4 - Passive - Neg - Up
%% column 5 - Passive - Pos - Down
%% column 6 - Passive - Pos - Up
%% column 7 - Passive - Overlap - &

%% column 1 - Attention - All - Down
%% column 2 - Attention - All - Up
%% column 3 - Attention - Neg - Down
%% column 4 - Attention - Neg - Up
%% column 5 - Attention - Pos - Down
%% column 6 - Attention - Pos - Up
%% column 7 - Attention - Overlap - &

%% Anatomical Lobes + Subcortical

%% line 1 - FRONTAL
%% line 48 - PARIETAL
%% line 63 - OCCIPITAL
%% line 77 - TEMPORAL
%% line 96 - SUBCORTICAL

line(1) = 1;
line(2) = 48;
line(3) = 63;
line(4) = 77;
line(5) = 96;
last_line = 118;

lobe{1} = 'FRONTAL';
lobe{2} = 'PARIETAL';
lobe{3} = 'OCCIPITAL';
lobe{4} = 'TEMPORAL';
lobe{5} = 'SUBCORTICAL';

nLobe = 5;

pcriterion = 0.01;

for iLobe=1:nLobe
    
    if iLobe < nLobe; start = line(iLobe); end_ = line(iLobe+1)-1; end;
    if iLobe == nLobe; start = line(iLobe); end_ = last_line; end;

    [H,P] = ttest(ResultsAnatomical(start:end_,1),ResultsAnatomical(start:end_,8));

    if P < pcriterion
        if mean(ResultsAnatomical(start:end_,1)) > mean(ResultsAnatomical(start:end_,8))
            disp(strcat(lobe{iLobe},':','Passive - All - Down - changes bigger'));
        else
            disp(strcat(lobe{iLobe},':','Attention - All - Down - changes bigger'));
        end
    end

    [H,P] = ttest(ResultsAnatomical(start:end_,2),ResultsAnatomical(start:end_,9));

    if P < pcriterion
        if mean(ResultsAnatomical(start:end_,2)) > mean(ResultsAnatomical(start:end_,9))
            disp(strcat(lobe{iLobe},':','Passive - All - Up - changes bigger'));
        else
            disp(strcat(lobe{iLobe},':','Attention - All - Up - changes bigger'));
        end
    end

    [H,P] = ttest(ResultsAnatomical(start:end_,3),ResultsAnatomical(start:end_,10));

    if P < pcriterion
        if mean(ResultsAnatomical(start:end_,3)) > mean(ResultsAnatomical(start:end_,10))
            disp(strcat(lobe{iLobe},':','Passive - Neg - Down - changes bigger'));
        else
            disp(strcat(lobe{iLobe},':','Attention - Neg - Down - changes bigger'));
        end
    end

    [H,P] = ttest(ResultsAnatomical(start:end_,4),ResultsAnatomical(start:end_,11));

    if P < pcriterion
        if mean(ResultsAnatomical(start:end_,4)) > mean(ResultsAnatomical(start:end_,11))
            disp(strcat(lobe{iLobe},':','Passive - Neg - Up - changes bigger'));
        else
            disp(strcat(lobe{iLobe},':','Attention - Neg - Up - changes bigger'));
        end
    end

    [H,P] = ttest(ResultsAnatomical(start:end_,5),ResultsAnatomical(start:end_,12));

    if P < pcriterion
        if mean(ResultsAnatomical(start:end_,5)) > mean(ResultsAnatomical(start:end_,12))
            disp(strcat(lobe{iLobe},':','Passive - Pos - Down - changes bigger'));
        else
            disp(strcat(lobe{iLobe},':','Attention - Pos - Down - changes bigger'));
        end
    end

    [H,P] = ttest(ResultsAnatomical(start:end_,6),ResultsAnatomical(start:end_,13));

    if P < pcriterion
        if mean(ResultsAnatomical(start:end_,6)) > mean(ResultsAnatomical(start:end_,13))
            disp(strcat(lobe{iLobe},':','Passive - Pos - Up - changes bigger'));
        else
            disp(strcat(lobe{iLobe},':','Attention - Pos - Up - changes bigger'));
        end
    end

    [H,P] = ttest(ResultsAnatomical(start:end_,7),ResultsAnatomical(start:end_,14));

    if P < pcriterion
        if mean(ResultsAnatomical(start:end_,7)) > mean(ResultsAnatomical(start:end_,14))
            disp(strcat(lobe{iLobe},':','Passive - Overlap - & - changes bigger'));
        else
            disp(strcat(lobe{iLobe},':','Attention - Overlap - & - changes bigger'));
        end
    end

end

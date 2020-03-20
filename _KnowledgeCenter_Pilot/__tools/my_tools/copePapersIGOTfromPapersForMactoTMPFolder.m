

%%% MOVE SELECTED PDFS TO TMP FOLDER

folder(1) = 'D:\__PDFs-Papers-for-Mac\PAPERS-FOR-MAC\PDF\';
folder(2) = 'D:\__PDFs-Papers-for-Mac\ENDNOTE\PDF\';
folder(3) = 'D:\__PDFs-Papers-for-Mac\READCUBE\PDF\';
folder(4) = 'D:\__PDFs-Papers-for-Mac\MISSING\PDF\';

tmp_folder = 'D:\__tmp\';

for iPMID=1:length(all_I_got)
    
    for iFolder=1:length(folder)

        system(sprintf('copy %s%s.pdf %s%s.pdf',folder{iFolder},int2str(all_I_got(iPMID)),tmp_folder,int2str(all_I_got(iPMID))));

    end

end
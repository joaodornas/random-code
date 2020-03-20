load sf_lsfv_dyn_forLab.mat
syms sf

warning off symbolic:sym:int:warnmsg1

for i=[4,6:12,14:ncells]
    tic
    for k=1:200
        disp(['computing LSFVs of cell #' num2str(i) ', time window ' num2str(k) ' (windowSize50)'])
        if ~isnan(pref_SFs_dyn_50(i,k))
            try
                fittedmodel=eval(['fittedmodels_50.cell' num2str(i) '.window' num2str(k)]);

                A=fittedmodel.A;
                bestSF=fittedmodel.ps;
                tuningWidth=fittedmodel.tuningWidth;
                skew=fittedmodel.skew;
                a=int((A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))))).*(log2(sf)/log2(16)-log2(bestSF)/log2(16)).^2,bestSF/16,bestSF);
                b=int((A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))))),bestSF/16,bestSF);

                try
                    LSFVs_dyn_50(i,k)=double(a)/double(b);

                catch

                    LSFVs_dyn_50(i,k)=nan;
                end

                clear a b A bestSF tuningWidth skew

            catch
                LSFVs_dyn_50(i,k)=nan;
            end
        end
    end

    clear fittedmodel
    save sf_lsfv_dyn_forLab
    disp(['elapsed time for cell ' num2str(i) ' = ' num2str(fixdec(toc/60,1)) ' minutes'])
end

for i=[6:12,14:ncells]
    tic
    for k=1:200
        disp(['computing LSFVs of cell #' num2str(i) ', time window ' num2str(k) ' (windowSize10)'])
        if ~isnan(pref_SFs_dyn_10(i,k))
            try
                fittedmodel=eval(['fittedmodels_10.cell' num2str(i) '.window' num2str(k)]);

                A=fittedmodel.A;
                bestSF=fittedmodel.ps;
                tuningWidth=fittedmodel.tuningWidth;
                skew=fittedmodel.skew;
                a=int((A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))))).*(log2(sf)/log2(16)-log2(bestSF)/log2(16)).^2,bestSF/16,bestSF);
                b=int((A*((exp((-((log2(sf)-log2(bestSF)).^2)./(2*(tuningWidth+skew.*(log2(sf)-log2(bestSF))).^2))))-(exp(-1/((skew)^2))))),bestSF/16,bestSF);

                try
                    LSFVs_dyn_10(i,k)=double(a)/double(b);

                catch
                    LSFVs_dyn_10(i,k)=nan;
                end
                clear a b A bestSF tuningWidth skew
            catch
                LSFVs_dyn_10(i,k)=nan;
            end
        end
    end

    clear fittedmodel
    save sf_lsfv_dyn_forLab
    disp(['elapsed time for cell ' num2str(i) ' = ' num2str(fixdec(toc/60,1)) ' minutes'])
end
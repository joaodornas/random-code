function plotFigureFitCorrelationEvent(Events_Conditions,Fits_Gauss_Conditions,name,plot_name,kind,folderpath,bar,Visible)

    condicoes = max(size(Fits_Gauss_Conditions,1),size(Fits_Gauss_Conditions,2));
    
    disp(name);
    
    for c=1:condicoes
        
        nEventos = max(size(Fits_Gauss_Conditions(c).eventos,1),size(Fits_Gauss_Conditions(c).eventos,2));

        if nEventos >=1

            for n=1:nEventos

                start = Events_Conditions(c).eventos(n).start;

                end_ = Events_Conditions(c).eventos(n).end;

                eval(strcat('auto_correlogram = Fits_Gauss_Conditions(',int2str(c),').eventos(',int2str(n),').correlogram_',kind,'.correlogram_norm;'));
                eval(strcat('auto_xcorr = Fits_Gauss_Conditions(',int2str(c),').eventos(',int2str(n),').xcorr_',kind,'.xcorr_norm;'));


                if mod(length(auto_xcorr),2) ~= 0

                        k = 1;
                        j = 0;

                else

                        k = 0;
                        j = 1;

                end

                if mod(length(auto_correlogram),2) ~= 0

                        w = 1;
                        u = 0;

                else

                        w = 0;
                        u = 1;


                end

                tracex = (-( length(auto_xcorr) - k )/2 + j):( ( length(auto_xcorr) - k )/2);
                trace = (- ( length(auto_correlogram) - w )/2 + u):( ( length(auto_correlogram) - w )/2);

                eval(strcat('A = Fits_Gauss_Conditions(',int2str(c),').eventos(',int2str(n),').correlogram_',kind,'.A;'))
                eval(strcat('Ax = Fits_Gauss_Conditions(',int2str(c),').eventos(',int2str(n),').xcorr_',kind,'.Ax;'));

                eval(strcat('s = Fits_Gauss_Conditions(',int2str(c),').eventos(',int2str(n),').correlogram_',kind,'.sigma;'))
                eval(strcat('sx = Fits_Gauss_Conditions(',int2str(c),').eventos(',int2str(n),').xcorr_',kind,'.sigmax;'));

                eval(strcat('mu = Fits_Gauss_Conditions(',int2str(c),').eventos(',int2str(n),').correlogram_',kind,'.mu;'))
                eval(strcat('mux = Fits_Gauss_Conditions(',int2str(c),').eventos(',int2str(n),').xcorr_',kind,'.mux;'));

                f = figure;

                x = trace;
                xx = tracex;

                y = A * exp((-(x - mu).^2)./(2*s^2));
                yy = Ax * exp((-(xx - mux).^2)./(2*sx^2));

                y = y ./ max(y);
                yy = yy ./ max(yy);

                set(gcf,'Visible',Visible);
                plot(trace,auto_correlogram,'r')
                hold on;
                set(gcf,'Visible',Visible);
                plot(x,y,'b');
                text(0.7*max(trace),0.9*max(auto_correlogram),strcat('SIGMA:',num2str(s)),'FontSize',20);
                text(0.7*max(trace),0.7*max(auto_correlogram),strcat('mu:',num2str(mu)),'FontSize',20);
                print(f,'-depsc',strcat(folderpath,name,bar,plot_name,'-','singleEventsAutoCorr-c',int2str(c),'-event-',int2str(n),'-sigma-',num2str(s),'-mu-',num2str(mu),'-START-',num2str(start),'-END-',num2str(end_),'-',kind,'.eps'));
                close all;

                g = figure;

                set(gcf,'Visible',Visible);
                plot(tracex,auto_xcorr,'r')
                hold on;
                set(gcf,'Visible',Visible);
                plot(xx,yy,'b');
                text(0.7*max(tracex),0.9*max(auto_xcorr),strcat('SIGMA:',num2str(sx)),'FontSize',20);
                text(0.7*max(tracex),0.7*max(auto_xcorr),strcat('mu:',num2str(mux)),'FontSize',20);
                print(g,'-depsc',strcat(folderpath,name,bar,plot_name,'-','singleEventsAutoCorr-c',int2str(c),'-event-',int2str(n),'-sigmax-',num2str(sx),'-mux-',num2str(mux),'-START-',num2str(start),'-END-',num2str(end_),'-',kind,'.eps'));
                close all;

            end

        end
        
    end


end


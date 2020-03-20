function [lat]=latency(psth,parameters,method)

%function [lat]=latency(psth,parameters,method)
%calculates latency from psth
%
%Inputs: psth is a structure containing psths for all stimulation
%conditions, which should be named condition1, condition2, ..., conditionn
%
%parameters is a structure containing the following fields:
%nconditions: number of stimulation conditions
%trial_duration: duration of trials in ms
%baseline_duration: duration of pre-stimulation period in ms
%stim_duration: duration of stimulation in ms
%analysis_period: period of evoked response to be used for latency
%calculation, in ms
%analysis_startTime: time after begining of stimulation after which
%analysis should start, in ms
%
%method is a string argument and should be set to 'derivative', for
%calculating latency from the maximal value of the psth's derivative;
%'gaussian_fit', for fitting the psth with the sum of 3 gaussians and
%settig latency as the peak of the first gaussian; 'percentile' for
%calculating latency as the 1st point at which psth exceeds the 95% percentile
%calculated from the baseline after 2 that have exceeded the 99% percentile; 
%'poisson' for a similar percentile method, which calculates it
%assuming a poisson distribution with lambda = mean(baseline); 'maxLikelihood' 
%for the maximum likelihood estimation method of Friedman & Priebe (1998) 
%J Neurosci Meth 83(2):185-194; and 'all' for using all 5 previous methods
%
%Output: lat is a vector containing latency for each condition, except 
%when all methods are used, in which case it will be a structure containing
%5 vectors, one for each method
%written by Pedro G. Vieira and Lucas Pinto in March 2009


%load(uigetfile ('*.mat'));

%get parameters
nconditions=parameters.nconditions;
baseline_duration=parameters.baseline_duration;
analysis_period=parameters.analysis_period;
analysis_startTime=parameters.analysis_startTime;

switch method

    case 'derivative'

        lat=zeros(nconditions,1);

        for i = 1:nconditions;

            psth_i=eval(['psth.condition' num2str(i) ';']);
            
            %calculate derivative and find peak
            derivative_psth = diff(psth_i(baseline_duration+1:(baseline_duration+analysis_period)));
            peak = find(derivative_psth == max(derivative_psth));

            lat(i) = peak(1);

        end

    case 'gaussian_fit'

        lat=zeros(nconditions,1);

        for i = 1:nconditions;

            psth_i=eval(['psth.condition' num2str(i) ';']);
            
            %fit psth with the sum of three gaussians
            psth_period=psth_i(baseline_duration:(baseline_duration+analysis_period));
            time=1:analysis_period+1;
            fresult = fit(time',psth_period','gauss3');
            
            %latency is the 1st peak
            lat(i) = fix(min([fresult.b1 fresult.b2 fresult.b3]));

        end

    case 'percentile'

        lat=zeros(nconditions,1);

        for i = 1:nconditions;

            psth_i=eval(['psth.condition' num2str(i)]);

            baseline=psth_i(1:baseline_duration);
            evoked_resp = psth_i(baseline_duration+analysis_startTime+1:baseline_duration+analysis_period);

            %calculate percentiles from experimental baseline distribution
            alpha1 = prctile(baseline,99);
            alpha2 = prctile(baseline,95);

            %find 1st point that exceeds 95% percentile after 2 that have
            %exceeded the 99% percentile
            for j=1:(size(evoked_resp,2)-2)
                if (evoked_resp(j) > alpha1 && evoked_resp(j+1)> alpha1 && evoked_resp(j+2) > alpha2)
                    lat(i)=j+analysis_startTime;
                    break
                else lat(i)=NaN;
                end
            end
        end

    case 'poisson'

        lat=zeros(nconditions,1);

        for i = 1:nconditions;

            psth_i=eval(['psth.condition' num2str(i)]);

            baseline=psth_i(1:baseline_duration).*1000;
            evoked_resp = psth_i(baseline_duration+analysis_startTime+1:baseline_duration+analysis_period).*1000;

            %calculate percentiles from poisson distribution with
            %lambda=mean(baseline)
            alpha1_poiss = poissinv(0.99,mean(baseline));
            alpha2_poiss = poissinv(0.95,mean(baseline));

            %find 1st point that exceeds 95% percentile after 2 that have
            %exceeded the 99% percentile
            for j=1:(size(evoked_resp,2)-2)
                if (evoked_resp(j) > alpha1_poiss && evoked_resp(j+1)> alpha1_poiss && evoked_resp(j+2) > alpha2_poiss)
                    lat(i)=j+analysis_startTime;
                    break
                else lat(i)=NaN;
                end
            end
        end
        
    case 'maxLikelihood'
        
        lat=zeros(nconditions,1);

        for i = 1:nconditions;

            psth_i=eval(['psth.condition' num2str(i)]);
            l=maxLikelihood_latency(psth_i,parameters);

            if isempty(l)==1
                lat(i)=NaN;
            else lat(i)=l;
            end

        end
        

    case 'all'

        lat.derivative=zeros(nconditions,1);
        lat.gaussian_fit=zeros(nconditions,1);
        lat.percentile=zeros(nconditions,1);
        lat.poisson=zeros(nconditions,1);
        lat.maxLikelihood=zeros(nconditions,1);

        for i = 1:nconditions;

            psth_i=eval(['psth.condition' num2str(i)]);

            %derivative
            derivative_psth = diff(psth_i(baseline_duration:(baseline_duration+analysis_period)));
            peak = find(derivative_psth == max(derivative_psth));

            lat.derivative(i) = peak;

            %gaussian
            psth_period=psth_i(baseline_duration:(baseline_duration+analysis_period));
            time=1:analysis_period+1;
            fresult = fit(time',psth_period','gauss3');

            lat.gaussian_fit(i) = fix(min([fresult.b1 fresult.b2 fresult.b3]));

            %percentile & poisson percentile
            baseline=psth_i(1:baseline_duration).*1000;
            evoked_resp = psth_i(baseline_duration+analysis_startTime+1:baseline_duration+analysis_period).*1000;

            alpha1 = prctile(baseline,99);
            alpha2 = prctile(baseline,95);
            alpha1_poiss = poissinv(0.99,mean(baseline));
            alpha2_poiss = poissinv(0.95,mean(baseline));

            for j=1:(size(evoked_resp,2)-2)
                if (evoked_resp(j) > alpha1 && evoked_resp(j+1)> alpha1 && evoked_resp(j+2) > alpha2)
                    lat.percentile(i)=j+analysis_startTime;
                    break
                else lat.percentile(i)=NaN;
                end
            end
            for j=1:(size(evoked_resp,2)-2)
                if (evoked_resp(j) > alpha1_poiss && evoked_resp(j+1)> alpha1_poiss && evoked_resp(j+2) > alpha2_poiss)
                    lat.poisson(i)=j+analysis_startTime;
                    break
                else lat.poisson(i)=NaN;
                end
            end
            
            %maximum Likelihood
            lat.maxLikelihood(i)=maxLikelihood_latency(psth_i,parameters);

        end
end
end

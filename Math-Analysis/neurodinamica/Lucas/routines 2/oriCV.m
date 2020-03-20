function [oriCV_sum,oriCV_int] = oriCV(responseVector,directionValues,fittedParameters,type)

%function: [oriCV_sum,oriCV_int]=oriCV(responseVector,directionValues,fittedParameters,type)
%cauculates orientation circular covariance using one of two methods or
%both: sum and integration, set in the type input: 'sum', 'integral' or
%'both'. responseVector and directionValues are needed for the sum method,
%and are two line vectors of the same size containing the neuron's
%responses to different directions and direction values in radians.
%fittedParamaters is a cftool object containing fitted model parameters
%(default matlab name fittedmodel1). 

switch type
    
    case 'sum'
        
        %orientation circular variance (with discrete sum)
        response=responseVector;
        oriCV_sum=1-(abs(sum(response.*exp(directionValues*2*i)))/(sum(response)));
        
        oriCV_int='not calculated';

    case 'integral'
        
        ndirections=size(directionValues,2);

        %model parameters
        baseline=fittedParameters.baseline;
        pref_A=fittedParameters.pref_A;
        pref_k=fittedParameters.pref_k;
        prefDir=fittedParameters.prefDir;
        anti_A=fittedParameters.anti_A;
        anti_k=fittedParameters.anti_k;
        antiPrefDir=fittedParameters.antiPrefDir;

        %orientation circular variance (with integration)
        syms x
        a=int((baseline+(pref_A*exp(pref_k*(cos(x-prefDir)-1)))+(anti_A*exp(anti_k*(cos(x-antiPrefDir)-1)))).*exp(x*2*i),x,directionValues(1),directionValues(ndirections));
        b=int(baseline+(pref_A*exp(pref_k*(cos(x-prefDir)-1)))+(anti_A*exp(anti_k*(cos(x-antiPrefDir)-1))),x,directionValues(1),directionValues(ndirections));
        oriCV_int=1-(double(abs(a))/double(b));
        
        oriCV_sum='not calculated';

    case 'both'
        
        ndirections=size(directionValues,2);
        
        %orientation circular variance (with discrete sum)
        response=responseVector;
        oriCV_sum=1-(abs(sum(response.*exp(directionValues*2*i)))/(sum(response)));
        
        %model parameters
        baseline=fittedParameters.baseline;
        pref_A=fittedParameters.pref_A;
        pref_k=fittedParameters.pref_k;
        prefDir=fittedParameters.prefDir;
        anti_A=fittedParameters.anti_A;
        anti_k=fittedParameters.anti_k;
        antiPrefDir=fittedParameters.antiPrefDir;
        
        %orientation circular variance (with integration)
        syms x
        a=int((baseline+(pref_A*exp(pref_k*(cos(x-prefDir)-1)))+(anti_A*exp(anti_k*(cos(x-antiPrefDir)-1)))).*exp(x*2*i),x,directionValues(1),directionValues(ndirections));
        b=int(baseline+(pref_A*exp(pref_k*(cos(x-prefDir)-1)))+(anti_A*exp(anti_k*(cos(x-antiPrefDir)-1))),x,directionValues(1),directionValues(ndirections));
        oriCV_int=1-(double(abs(a))/double(b));
end

end
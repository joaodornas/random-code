function [histogram_values,histogram_Xaxis]=makeHistogram(datavector,binwidth,plot)

%[histogram_values,histogram_Xaxis]=makeHistogram(datavector,binwidth,plot)
%
%returns histogram values and axis for plotting with the 'bar' plot
%command. Inputs are data vector, which can be either line or column, and
%desired binwidth of the histogram. Bin limits are automatically rounded to 
%nearest value in the binwidth scale. Plot is a string command and should 
%be set to 'on' (default) to automatically plot the histogram.
%
%by Lucas Pinto in April 2009

if nargin<2
    'error! I need at least 2 inputs: datavector and binwidth'
elseif nargin==2
    plot='on';
end

%find minimum and maximum values in data vector
minvalue=min(datavector);
maxvalue=max(datavector);

%calculate bins
histogram_bins=round(minvalue)-binwidth:binwidth:round(maxvalue)+binwidth;
nbins=numel(histogram_bins);
histogram_values=zeros(1,nbins);

%fill histogram
for i=1:nbins
    a=find(datavector >= histogram_bins(i) & datavector < histogram_bins(i)+binwidth);
    histogram_values(i)=numel(a);
end

%axis for plotting with the 'bar' command
histogram_Xaxis=histogram_bins+(binwidth/2);

%plot histogram
switch plot
    case 'on'
        figure;
        bar(histogram_Xaxis,histogram_values);
end
end

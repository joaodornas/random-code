function geom_std=geostd(datavector)

%calculates geometric standard deviation of a 1D vector, ignoring nan's 

a=isempty(isnan(datavector));

if a==0
    datavector=datavector(~isnan(datavector));
end


geom_std=exp(sqrt(sum((log(datavector)-log(geomean(datavector))).^2)/(numel(datavector)-1)));
function zscored = nanzscore(timeseries)

    zscored = ( timeseries - mean(timeseries(~isnan(timeseries))) ) ./ std(timeseries(~isnan(timeseries)));

end
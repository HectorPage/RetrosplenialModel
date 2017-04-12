function tuning_curve = generate_tuning_curve(trueHD,Rates,nBins,binEdges)
%generate_tuning_curve -  Generates a tuning curve for a given cell, from
%model data

%Created by and copyright held by Dr. Hector JI Page 24/10/16

%% Binning firing data and taking mean
[~,whichBin] = histc(trueHD,binEdges);

binMean = zeros(nBins,1);
for index = 1:nBins
    flagBinMembers = (whichBin == index);
    binMembers = Rates(flagBinMembers);
    binMean(index) = nanmean(binMembers);
    if(isnan(binMean(index)))
        binMean(index) = 0;
    end
end

%% SMOOTHING ACROSS 0-360 divide

tuning_curve = zeros(numel(binMean)+1,1);

%Doubling the final element
tuning_curve(1:end-1) = binMean;
tuning_curve(end) = binMean(1);
if sum(tuning_curve)>0.01
    %Smoothing tuning curve
    tuning_curve = smooth(tuning_curve)';
else
    tuning_curve = tuning_curve';
end
%Normalising tuning curve
%tuning_curve = tuning_curve/max(tuning_curve)';


end


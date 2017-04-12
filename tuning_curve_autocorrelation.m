function [autocorrelation, pks, separations] = tuning_curve_autocorrelation(threshold,tuning_curve)
%tuning_curve_autocorrelation - autocorrelates a tuning curve to detect
%planes of symmetry and numbers of peaks in tuning curve. Used for bipolar
%and tripolar cells

%Counts peaks in the autocorrelations, and the locations of these peaks in
%order to calculate their separation

%Created by and copyright held by Dr. Hector JI Page 24/10/16

autocorrelation = cconv(tuning_curve,conj(fliplr(tuning_curve)),numel(tuning_curve));
separations = zeros(2,1);

[pks, locs] = findpeaks(autocorrelation);
pks = pks(pks>threshold);
locs = locs(pks>threshold);
pks = numel(pks) +1; %num_peaks +1, as doesn't detect peak at 0
if numel(pks)>1
    separations(1) = abs(locs(2) - locs(1)) * 6;
    if(numel(locs)>2)
        separations(2) = abs(locs(3)-locs(2)) * 6;
    end
end


end


function PVector = calculate_PVector(Rates,PFDs)

%Calculates the population vector of a set of cells with firing rate =
%rates and preferred firing directions = PFDs.

%Created by and copyright held by Dr. Hector JI Page 03/11/16

fav_view_sin = sind(PFDs)'; %sine
fav_view_cos = cosd(PFDs)'; %cosine

vector_1 = sum(bsxfun(@times,Rates,fav_view_sin),1);
vector_2 = sum(bsxfun(@times,Rates,fav_view_cos),1);

PVector = atan2d(vector_1,vector_2);

PVector(PVector<=0) = PVector(PVector<=0) + 360;

end
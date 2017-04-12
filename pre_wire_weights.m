function [weight_profile, unconnected] = pre_wire_weights(model_region,num_presynaptic,num_postsynaptic,pfd_pre,pfd_post,sigma)
%This function pre-wires a Gaussian weight profile between two sets of
%neurons, based on similarity of PFDs.

%NOTE - THIS FUNCTION WILL NOT PRE-WIRE BIPOLAR OR TRIPOLAR WEIGHT PROFILES

%Created by and copyright held by Dr. Hector JI Page 26/09/16

if(strcmpi(model_region,'HD_to_conjunctive'))
    
    %For HD to conjunctive, the latter half of conjunctive cells do not
    %have any weight
    
    pfd_pre = pfd_pre';    
    pfd_pre = pfd_pre(:,ones(num_postsynaptic,1),:);
    pfd_post = pfd_post(ones(num_presynaptic,1),:);
    
    diff = abs(pfd_pre-pfd_post);
    diff(diff>180) = 360 - diff(diff>180);
    diff = diff./sigma;
    diff = diff.^2;
    diff = diff.*-0.5;
    diff = exp(diff);
    weight_profile = diff .*0.5;
    
    %Now set last half of postsynaptic weights to be zero
    weight_profile(:,num_presynaptic+1:end) = 0;
    unconnected = num_presynaptic+1:num_postsynaptic; % not connected
    
 else
    %For the rest of the areas connectivity wired based on PFDs (where
    %conj PFDs are set up to be bimodal prior to calling this function)!
    
    pfd_pre = pfd_pre';
    pfd_pre = pfd_pre(:,ones(num_postsynaptic,1),:);
    pfd_post = pfd_post(ones(num_presynaptic,1),:);
    
    diff = abs(pfd_pre-pfd_post);
    diff(diff>180) = 360 - diff(diff>180);
    diff = diff./sigma;
    diff = diff.^2;
    diff = diff.*-0.5;
    diff = exp(diff);
    weight_profile = diff .*0.5;
    
    unconnected = [];
end


end
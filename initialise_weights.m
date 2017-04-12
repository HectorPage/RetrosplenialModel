function [weight_profile,unconnected] = initialise_weights(model_region,num_presynaptic,num_postsynaptic,type,proportion)
% Generates a random weight profile
% type = 'incomplete' means proportion (in range [0,1]) postsynaptic don't receive from
% presynaptic 
% type = 'complete' means full interconnectivity

%proportion is in % e.g. 50

if(strcmpi(type,'complete'))

    weight_profile = rand(num_presynaptic,num_postsynaptic);
    unconnected = [];
    
elseif(strcmpi(type,'incomplete'))
    
    
    proportion = proportion/100;
    
     if(strcmpi(model_region,'HD_to_conjunctive'))
         weight_profile = rand(num_presynaptic,num_postsynaptic);
         
          %Now set last half of postsynaptic weights to be zero
          weight_profile(:,num_presynaptic+1:end) = 0;
          unconnected = num_presynaptic+1:num_postsynaptic; % not connected
         
     else
         weight_profile = rand(num_presynaptic,num_postsynaptic);
         unconnected = datasample(1:num_postsynaptic,(num_postsynaptic*proportion),'Replace',false);
         weight_profile(:,unconnected) = 0; %set unconnected weights to 0 - will have to edit code to keep them at 0
     end
   
    
else
   disp('Error: Usage of type must be ''incomplete'' or ''complete'''); 
end

weight_profile = weight_profile./10; %make initial weights pretty small

end


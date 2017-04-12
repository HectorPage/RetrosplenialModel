function interrogate_vis_input(filename)
% This function tests a full SO model after simulation, to see how the vis
%landmarks can drive specific head directions (or not)

%Created by and copyright held by Dr. Hector JI Page 18/01/17

disp('Now testing landmark influence');

%% READ PARAMETERS
param_path = strcat(pwd,'/',filename);

[~, ~, hd_cells,ADN_cells,conjunctive_cells, num_landmarks, phi_HDconjunctive, phi_conjunctiveHD, phi_vis,...
    phi_ADNHD,phi_HDADN,phi_ADNADN,tau_HD, tau_conjunctive, tau_ADN, timestep_size, HD_inhibition, conjunctive_inhibition, ADN_inhibition,...
    conduction_delay,~, ~, ~, ~, ~,~, ~, ~,~,~,~,~,~] = read_RSC_v2_params(param_path);

time = num_landmarks; % ten seconds of testing

%% Generate model terms needed
epochs = time/timestep_size;
HD_coefficient = timestep_size/tau_HD;
conjunctive_coefficient = timestep_size/tau_conjunctive;
ADN_coefficient = timestep_size/tau_ADN;

HD_conjunctive_scale = phi_HDconjunctive/hd_cells;
conjunctive_HD_scale = phi_conjunctiveHD/conjunctive_cells;
vis_conjunctive_scale = phi_vis/num_landmarks;
ADN_HD_scale = phi_ADNHD/ADN_cells;
HD_ADN_scale = phi_HDADN/hd_cells;
ADN_ADN_scale = phi_ADNADN/ADN_cells;

disp('Model terms computed');

%% Read in weights
fid = fopen('vis_to_conjunctive_weights_final.bdat', 'rb');
vis_to_conjunctive_weights = fread(fid, [num_landmarks,conjunctive_cells], 'float32')';
vis_to_conjunctive_weights = vis_to_conjunctive_weights';
fclose(fid);

fid = fopen('HD_to_conjunctive_weights_final.bdat', 'rb');
HD_to_conjunctive_weights = fread(fid, [hd_cells,conjunctive_cells], 'float32')';
HD_to_conjunctive_weights = HD_to_conjunctive_weights';
fclose(fid);

fid = fopen('conjunctive_to_HD_weights_final.bdat', 'rb');
conjunctive_to_HD_weights = fread(fid, [conjunctive_cells,hd_cells], 'float32')';
conjunctive_to_HD_weights = conjunctive_to_HD_weights';
fclose(fid);

fid = fopen('ADN_to_HD_weights.bdat', 'rb');
ADN_to_HD_weights = fread(fid,[hd_cells,hd_cells],'float32');
ADN_to_HD_weights = ADN_to_HD_weights';
fclose(fid);

fid = fopen('HD_to_ADN_weights.bdat', 'rb');
HD_to_ADN_weights = fread(fid,[hd_cells,hd_cells],'float32');
HD_to_ADN_weights = HD_to_ADN_weights';
fclose(fid);

fid = fopen('ADN_to_ADN_weights.bdat', 'rb');
ADN_to_ADN_weights = fread(fid,[hd_cells,hd_cells],'float32');
ADN_to_ADN_weights = ADN_to_ADN_weights';
fclose(fid);


%% Allocate memory to store rates etc.
prev_rates_HD = zeros(hd_cells,1);
prev_activations_HD = zeros(hd_cells,1);

prev_rates_conjunctive = zeros(conjunctive_cells,1);
prev_activations_conjunctive = zeros(conjunctive_cells,1);

prev_rates_ADN = zeros(ADN_cells,1);
prev_activations_ADN = zeros(ADN_cells,1);

%Here's the rates over time - saving every 100th timestep to reduce memory
%usage

rates_HD_time = zeros(hd_cells,epochs/100);
rates_conjunctive_time = zeros(conjunctive_cells,epochs/100);
rates_vis_time = zeros(num_landmarks,epochs/100);

rates_ADN_time = zeros(ADN_cells,epochs/100);

disp('Network memory allocated');


%% Loading conduction buffers
BufferLength = conduction_delay/timestep_size;
circBuff_HD = file_load(BufferLength, hd_cells, 'HDBuffer.bdat');
circBuff_conjunctive = file_load(BufferLength, conjunctive_cells, 'CONJBuffer.bdat');
circBuff_vis =  file_load(BufferLength, num_landmarks, 'VISBuffer.bdat');
circBuff_ADN = file_load(BufferLength, ADN_cells, 'ADNBuffer.bdat');

%Put them the right way around
circBuff_HD = circBuff_HD';
circBuff_conjunctive = circBuff_conjunctive';
circBuff_vis = circBuff_vis';
circBuff_ADN = circBuff_ADN';

current_landmark = 1;
one_sec = 1/timestep_size;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN MODEL TEST
for timestep = 1:epochs
    
    disp(['Timestep: ',num2str(timestep)]);
    
    %% NO PI input
    PI_input = zeros(hd_cells,1); %no PI input to ADN, just seeing how RSC HD influences it...
    
    %% Set landmark cell firing
    if (~mod(timestep,one_sec)) && timestep< one_sec*num_landmarks  %change landmark every second
        current_landmark = current_landmark +1;
    end
    
    rates_vis = zeros(num_landmarks,1);
    rates_vis(current_landmark) = 1;
    
    %% Calculate activations
    
    %ADN
    ADN_inhibitory = sum(prev_rates_ADN*ADN_inhibition);
    
    HD_input_to_ADN = sum(bsxfun(@times,HD_to_ADN_weights,circBuff_HD(:,end)),1)';
    
    ADN_input_to_ADN = sum(bsxfun(@times,ADN_to_ADN_weights,circBuff_ADN(:,end)),1)';
    
    activations_ADN = (1.0 - ADN_coefficient).*prev_activations_ADN...
        + ADN_coefficient.*(HD_input_to_ADN.*HD_ADN_scale)...
        + ADN_coefficient.*(ADN_input_to_ADN.*ADN_ADN_scale)...
        + ADN_coefficient.*PI_input...
        - ADN_coefficient.*(ADN_inhibitory.*(1/ADN_cells));
    
    %HD
    HD_inhibitory = sum(prev_rates_HD.*HD_inhibition);
    
    %Potentially faster bsxfun version
    conjunctive_input = sum(bsxfun(@times,conjunctive_to_HD_weights,circBuff_conjunctive(:,end)),1)';
    
    ADN_input = sum(bsxfun(@times,ADN_to_HD_weights,circBuff_ADN(:,end)),1)';
    
    %ADN_input = zeros(hd_cells,1);
    activations_HD = (1.0 - HD_coefficient).*prev_activations_HD...
        + HD_coefficient.*(ADN_input.*ADN_HD_scale)...
        + HD_coefficient.*(conjunctive_input.*conjunctive_HD_scale)...
        - HD_coefficient.*(HD_inhibitory.*(1/hd_cells));
    
    
    %Conjunctive
    conjunctive_inhibitory = sum(prev_rates_conjunctive.*conjunctive_inhibition);
    
    %Potentially faster bsxfun version
    HD_input = sum(bsxfun(@times,HD_to_conjunctive_weights,circBuff_HD(:,end)),1)';
    %HD_input = zeros(conjunctive_cells,1);
    
    %Potentially faster bsxfun version
    vis_input = sum(bsxfun(@times,vis_to_conjunctive_weights,circBuff_vis(:,end)),1)';
    
    activations_conjunctive = (1.0 - conjunctive_coefficient).*prev_activations_conjunctive...
        + conjunctive_coefficient.*(vis_input.*vis_conjunctive_scale)...
        + conjunctive_coefficient.*(HD_input.*HD_conjunctive_scale)...
        - conjunctive_coefficient.*(conjunctive_inhibitory.*(1/conjunctive_cells));
    
    %% Calculate rates
    
    %ADN
    rates_ADN = tanh(activations_ADN);
    rates_ADN(rates_ADN<0.0) = 0.0;
    
    %HD
    rates_HD = tanh(activations_HD); %tanh transfer function - can change if needed
    rates_HD(rates_HD<0.0) = 0.0;
    
    %Conjunctive
    rates_conjunctive = tanh(activations_conjunctive); %tanh transfer function - can change if needed
    rates_conjunctive(rates_conjunctive<0.0) = 0.0;
    
    %% Record previous states
    %Rates and acts
    prev_rates_ADN = rates_ADN;
    prev_activations_ADN = activations_ADN;
    
    prev_rates_HD = rates_HD;
    prev_activations_HD = activations_HD;
    
    prev_rates_conjunctive = rates_conjunctive;
    prev_activations_conjunctive = activations_conjunctive;
    
    %prev_rates_vis = rates_vis; %unused with buffer, but will be used in simulations without buffers
    
    %Rates and acts over time
    if(~mod(timestep,100)) %i.e. if timestep divisible by 100
        rates_HD_time(:,timestep/100) = rates_HD;
        rates_conjunctive_time(:,timestep/100) = rates_conjunctive;
        rates_vis_time(:,timestep/100) = rates_vis;
        
        rates_ADN_time(:,timestep/100) = rates_ADN;
        
    end
    
    %% Fill buffers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Buffer is set up so first element is most recent  %
    % and last timestep is the rates delta t ago, where %
    %     delta t is the conduction delay used....      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    circBuff_HD = [rates_HD circBuff_HD(:,1:end-1)]; %drops end value, and loads most recent value in at start
    circBuff_conjunctive = [rates_conjunctive circBuff_conjunctive(:,1:end-1)];
    circBuff_vis = [rates_vis circBuff_vis(:,1:end-1)];
    circBuff_ADN = [rates_ADN circBuff_ADN(:,1:end-1)];
    
    
    
end
disp('Simulation finished');


%% Save test rates
fileID = fopen('visRates_interrogation.bdat','wb');
fwrite(fileID,rates_vis_time,'float32');
fclose(fileID);

fileID = fopen('HDRates_interrogation.bdat','wb');
fwrite(fileID,rates_HD_time,'float32');
fclose(fileID);

fileID = fopen('conjunctiveRates_interrogation.bdat','wb');
fwrite(fileID,rates_conjunctive_time,'float32');
fclose(fileID);

fileID = fopen('ADNRates_interrogation.bdat','wb');
fwrite(fileID,rates_ADN_time,'float32');
fclose(fileID);

disp('Firing Rate Data Saved');

pause;
%% Plot test rates
interrogation_movie(time,timestep_size,conjunctive_cells,hd_cells,hd_cells,num_landmarks,0, time)


end

function rates = file_load(cells, steps, fname)

rates = zeros(cells, steps);

fid = fopen(fname, 'rb');

rates = fread(fid, [steps, cells], 'float32')';

fclose(fid);


end
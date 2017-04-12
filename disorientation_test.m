function disorientation_test(filename,initial_position,initial_heading,visual_scene)
% This function tests a full SO model after simulation, to see how the vis
% to CONJ to HD weights allow the system to correct HD with
% disorientated animal (i.e. ADN input at start at random value)

%visual_scene is an option for the 12 landmark, 6*12 visual cell version of
%a full visual scene influencing the animal

%it is updated to have a full ADN network, v1 did not have this

%Created by and copyright held by Dr. Hector JI Page 03/02/17


disp('Performing disorientation test');


%% READ PARAMETERS
param_path = strcat(pwd,'/',filename);

[~, ~, hd_cells,ADN_cells,conjunctive_cells, num_landmarks, phi_HDconjunctive, phi_conjunctiveHD, phi_vis,...
    phi_ADNHD,phi_HDADN,phi_ADNADN,tau_HD, tau_conjunctive, tau_ADN, timestep_size, HD_inhibition, conjunctive_inhibition, ADN_inhibition,...
    conduction_delay,~, ~, ~, ~, ~,...
    ~, ~, ~,PI_noise_percentage,PI_noise_width, PI_strength, PI_sigma, ~] = read_RSC_v2_params(param_path);

time = 30.0; % thirty seconds of testing

if visual_scene
    num_landmarks = 12*6;
end

% phi_conjunctiveHD = phi_conjunctiveHD * 6;
% phi_ADNHD = phi_ADNHD * 0.75;
% phi_HDADN = phi_HDADN * 2;


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

disp('Model terms computed')

%% Setting up landmark coordinates

if visual_scene %putting the clock face landmarks into the simulation
    num_landmarks = 12; %returning to 12 to set up the 12 clock number coordinates
    
    disp('Landmarks set to 12 clock numbers')
    
    increment = 2*pi/num_landmarks;
    angles = increment:increment:2*pi;
    initial_coords =[cos(angles);sin(angles)];
    initial_coords = initial_coords./2; %0.5m radius circle, centered on 0,0
    
    landmark_coords = zeros(2,12);
    for idx = 1:12 %reorder these  to be a clock face
        newindex = idx + 2;
        newindex(newindex>12) = newindex(newindex>12) - 12;
        landmark_coords(:,idx) = initial_coords(:,newindex);
    end
    
    landmark_coords = fliplr(landmark_coords); %now finally index corresponds with clock face position
    
    landmark_coords = landmark_coords'; %same dimensions as before
    
    num_landmarks = 12*6; %there are a max of 6 landmarks seen at a time with a 90 deg field of view
    %landmarks must be seen in order e.g. can't see 1,8,and 4
    %therefore the 6 combos are n, (n, n+1), (n, n+1, n+2) etc up to n+6.
    
    LandmarksVisible = zeros(12,epochs/100);
    
else
    if num_landmarks == 12
        %landmark_coords = 1.2*rand(num_landmarks,2); %x rows, 2 columns
        %SETTING LANDMARKS TO BE AROUND EXTERIOR OF COMPARTMENT, NON-RANDOM
        landmark_coords = [0 ,0;
            0.4,0;
            0.8,0;
            1.2,0;
            0,0.4;
            0,0.8;
            0,1.2;
            0.4,1.2;
            0.8,1.2;
            1.2,1.2;
            1.2,0.8
            1.2,0.4];
    elseif num_landmarks == 2 %these 2 landmarks are the edges of the cue card in the circular arena
        landmark_coords = [0.4,1.2;
            0.8, 1.2];
    end
end

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

%% Allocating memory for normalisation operations
headings_vector = zeros((epochs/100)+1,1);
headings_matrix = zeros((epochs/100)+1,2);

disp('Network memory allocated');


%% Do a short random walk
if visual_scene
    initial_position = [0,0]; %different as origin of circle
    sigma_head_turns = 15; %this chosen as works better for the circular arena simulations
    step_scale_lower = 0.05;
    step_scale_upper = 0.05;
    [headings,positions] = circular_arena_walk_v3('disorientation',epochs/100, initial_position,initial_heading,sigma_head_turns, step_scale_lower,step_scale_upper);
    disp('Circular arena walk finished');
else
    sigma_head_turns = 25;
    step_scale_lower = 0.05;
    step_scale_upper = 0.05;
    [headings,positions] = one_compartment_walk_v2('disorientation',1,epochs/100, initial_position,initial_heading,sigma_head_turns, step_scale_lower,step_scale_upper);
    disp('One compartment walk finished');
end
%% Get true HD (in deg) over course of simulation
[headings, ~, ~] = normalise_matrix(headings,'r',headings_vector,headings_matrix); %normalise each row of headings
trueHD = atan2(headings(:,2),headings(:,1));
trueHD(trueHD<0) = trueHD(trueHD<0) + (2*pi);
trueHD = trueHD * (180/pi); %convert to degrees


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


%% Setting HD PFDS
[HD_PFDs] = set_PFDs(hd_cells);

%% Converting percentages

    PI_noise_percentage = PI_noise_percentage/100;
    PI_noise_width = PI_noise_width/100;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN MODEL TEST
for timestep = 1:epochs
    
    disp(['Timestep: ',num2str(timestep)]);
    
    %% Set PI input
    if(timestep<2)
        PI_direction = trueHD(timestep); %first timestep ADN activity is initialised to true HD
        
        %now add random disorientation to starting HD...
        PI_direction = PI_direction + (360*rand(1,1));
        PI_direction(PI_direction>360) = PI_direction(PI_direction>360) - 360; %keep it bounded
        
    else
        %All other timesteps, ADN path integration input is noisy
        if trueHD(ceil(timestep/100))~=trueHD(ceil((timestep-1)/100)); %only a PI increment when trueHD changes
            %Initial increment is change in true HD
            PI_increment = atan2d(sind(trueHD(ceil(timestep/100))-trueHD(ceil((timestep-1)/100))),...
                cosd(trueHD(ceil(timestep/100))-trueHD(ceil((timestep-1)/100))));
            
            %Add noise, drawn from Guassian probility cetnered on increment
            %minus a percentage of increment, with sigma/width a certain
            %percentage of increment
            
            mean = PI_increment - (PI_noise_percentage*PI_increment);
            sigma = PI_noise_width*abs(PI_increment);
            
            PI_increment = normrnd(mean,sigma);
            
            %Calculate ADN PVector
            
            ADN_PV = calculate_PVector(prev_rates_ADN, HD_PFDs);
            
            
            %New PI input direction is noisy path integration increment from
            %ADNPI
            PI_direction = ADN_PV + PI_increment;
            
            PI_direction(PI_direction>360) = PI_direction(PI_direction>360) - 360;
            PI_direction(PI_direction<=0) = PI_direction(PI_direction<=0) + 360;
        end
    end
    
    distance1 = abs((ones(1,hd_cells)*PI_direction) - HD_PFDs); %heading and position and true HD change every 100 timesteps....
    distance2 = (ones(1,hd_cells)*360.0) - distance1;
    distance = min([distance1;distance2],[],1); %wraparound distance
    distance = (distance./PI_sigma).^2;
    
    %Guassian input, based on similarity of trueHD to PFD of each cell
    PI_input = PI_strength .* exp(-0.5.*distance);
    PI_input = PI_input'; %making correct dimension for later usage
 
    
    
    %% Set landmark cell firing
    %landmark vectors
    landmark_vectors = landmark_coords - repmat(positions(ceil(timestep/100),:),size(landmark_coords,1),1);
    
    %angle between landmark vector and HD
    landmark_vectors = normalise_matrix(landmark_vectors,'r');
    
    curr_heading = headings(ceil(timestep/100),:);
    if visual_scene
        dot_prods = dot(landmark_vectors, curr_heading(ones(12,1),:),2);
    else
        dot_prods = dot(landmark_vectors, curr_heading(ones(num_landmarks,1),:),2);
    end
    %I think the landmark angles need to use the acos(dotproduct)
    %method here
    landmark_angles = acos(dot_prods); %angle between vectors a and b = acos(dot(a,b));
    
    if(~visual_scene)
        %if angle <45, landmark within 45 deg field of view
        rates_vis = landmark_angles < pi/4; %(N.B. radians) %vis cells binary landmark detectors
    else
        
        visible_landmarks = landmark_angles < pi/4;
        
        %Here working out which combos of landmarks and setting the
        %correct cells
        
        %1. get indices of visible landmarks
        indices = find(visible_landmarks);
        
        if numel(indices)>6
            error('Landmark visibility error');
        end
        if isempty(indices)
            rates_vis = zeros(num_landmarks,1);
        else
            %2. find index where landmark pattern starts %(i.e. clockface order with no discontinuity)
            if numel(indices)<2
                start = indices(1);
            else
                
                if numel(indices)< sum(diff(indices))+1; %if differences between adjacents aren't all 1
                    start = indices(find(diff(indices)>1) + 1); %start of combined landmark is
                else
                    start = indices(1);
                end
            end
            
            %3.vis cell index = initial landmark - 1 * 6 +
            %numel(indices)
            vis_cell_index = ((start-1) *6) + numel(indices);
            
            %NOW USE TO DETERMINE RATES
            rates_vis = zeros(num_landmarks,1);
            rates_vis(vis_cell_index) = 1; %set the one visual cell that represents this landmark combination to fire
        end
    end
    
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
    
    activations_HD = (1.0 - HD_coefficient).*prev_activations_HD...
        + HD_coefficient.*(ADN_input.*ADN_HD_scale)...
        + HD_coefficient.*(conjunctive_input.*conjunctive_HD_scale)...
        - HD_coefficient.*(HD_inhibitory.*(1/hd_cells));
    
    
    %Conjunctive
    conjunctive_inhibitory = sum(prev_rates_conjunctive.*conjunctive_inhibition);
    
    
    %Potentially faster bsxfun version
    HD_input = sum(bsxfun(@times,HD_to_conjunctive_weights,circBuff_HD(:,end)),1)';
    
    
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
    prev_activations_conjunctive =activations_conjunctive;
    
    %prev_rates_vis = rates_vis; %unused with buffer, but will be used in simulations without buffers
    
    %Rates and acts over time
    if(~mod(timestep,100)) %i.e. if timestep divisible by 100
        rates_HD_time(:,timestep/100) = rates_HD;
        rates_conjunctive_time(:,timestep/100) = rates_conjunctive;
        rates_vis_time(:,timestep/100) = rates_vis;
        
        rates_ADN_time(:,timestep/100) = rates_ADN;
        
        if visual_scene
            LandmarksVisible(:,timestep/100) = visible_landmarks;
        end
        
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
fileID = fopen('visRates_disor.bdat','wb');
fwrite(fileID,rates_vis_time,'float32');
fclose(fileID);

fileID = fopen('HDRates_disor.bdat','wb');
fwrite(fileID,rates_HD_time,'float32');
fclose(fileID);

fileID = fopen('conjunctiveRates_disor.bdat','wb');
fwrite(fileID,rates_conjunctive_time,'float32');
fclose(fileID);

fileID = fopen('ADNRates_disor.bdat','wb');
fwrite(fileID,rates_ADN_time,'float32');
fclose(fileID);

%Save landmarks
%save landmark visibility if needed
if visual_scene
    fileID = fopen('LandmarkVisibility_disor.bdat','wb');
    fwrite(fileID,LandmarksVisible,'float32');
    fclose(fileID);
end


disp('Firing Rate Data Saved');

%% Plot test rates
RSC_model_plot_v3('disorientation','simulation',1,time, timestep_size, conjunctive_cells,num_landmarks, hd_cells, 'save',...
    rates_vis_time, rates_HD_time, rates_conjunctive_time, rates_ADN_time);

end

function rates = file_load(cells, steps, fname)

rates = zeros(cells, steps);

fid = fopen(fname, 'rb');

rates = fread(fid, [steps, cells], 'float32')';

fclose(fid);


end
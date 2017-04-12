function simulate_RSC_v3(model_type,filename)
%RSC MODEL RUN WITh VARYING MODEL TYPE:

%THIS IS VERSION 3: includes a full ADN network, with a kalman filtering of
%true HD - ADN input is a gaussian with location 'previous ADN location plus change in
%true HD plus noise', and has added to it a Gaussian input from RSC

%Above is version 2 edits, version 3 now has 12 clock-face landmarks spaced
%around the exterior of the environment, comibnations of landmarks are
%processed differently

%'full SO' = no pre-wired connectivity with learning of landmarks, one
%compartment random walk model_mode = 1;

%'cue conflict' = pre-wired vision, but the visual relationship to HD
%changes to generate conflict, one compartment random walk model_mode = 2;

%'bipolar' = pre-wired vision, HD to conj and conj to HD self-organises,
%switches based on compartment with two-compartment walk model_mode = 3;

%'tripolar' = pre-wired vision, HD to conj and conj to HD self-organises,
%switches based on compartment with three-compartment walk model_mode = 4;

%'filename' = the name of the paramns excel file, endining in .xlsx
%e.g. 'params.xlsx'


%Created by and copyright held by Dr. Hector JI Page 22/09/16


close all; %get rid of figures

model_type = model_type(1); %first letter is different, so just use that for speed

if(~strcmpi(model_type,'b')&& ~strcmpi(model_type,'t')...
        &&~strcmpi(model_type,'c')&&~strcmpi(model_type,'f')&&~strcmpi(model_type,'v'))
    error('Usage must be one of: \n ''full SO'' \n \n ''bipolar'' \n \n ''tripolar'' \n \n ''cue conflict'' \n\ ''visual scene''');
end

%% READ PARAMETERS
param_path = strcat(pwd,'/',filename);

[time, training_time, HD_cells,ADN_cells,conjunctive_cells, num_landmarks, phi_HDconjunctive, phi_conjunctiveHD, phi_vis,...
    phi_ADNHD,phi_HDADN,phi_ADNADN,tau_HD, tau_conjunctive, tau_ADN, timestep_size, HD_inhibition, conjunctive_inhibition, ADN_inhibition,...
    conduction_delay,learning_rate, vis_conjunctive_weight_sigma, HD_conjunctive_weight_sigma, conjunctive_HD_weight_sigma, vis_sigma,...
    ADN_ADN_weight_sigma, ADN_HD_weight_sigma, HD_ADN_weight_sigma,PI_undershoot_percentage,PI_noise_width, PI_strength, PI_sigma, LTD] = read_RSC_v2_params(param_path);

if(~strcmpi(model_type,'f')) %if not fully self-organising in terms of vis to conj
    num_landmarks = HD_cells;
end

if(strcmpi(model_type,'c'))
    conflict_size = 60;
end

if(strcmpi(model_type,'v'))
    num_landmarks = 12*6;
end


disp(['Paramaters read in from: ', param_path]);

%% Generate model terms needed
epochs = time/timestep_size;
HD_coefficient = timestep_size/tau_HD;
conjunctive_coefficient = timestep_size/tau_conjunctive;
ADN_coefficient = timestep_size/tau_ADN;

HD_conjunctive_scale = phi_HDconjunctive/HD_cells;
conjunctive_HD_scale = phi_conjunctiveHD/conjunctive_cells;
vis_conjunctive_scale = phi_vis/num_landmarks;
ADN_HD_scale = phi_ADNHD/ADN_cells;
HD_ADN_scale = phi_HDADN/HD_cells;
ADN_ADN_scale = phi_ADNADN/ADN_cells;

disp('Model terms computed')

%% Set up landmarks
if(strcmpi(model_type,'v'))
    
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
    
    formatSpec = '%f\n';
    
    fileID = fopen('landmarks_x.txt','w');
    fprintf(fileID, formatSpec, landmark_coords(:,1));
    fclose(fileID);
    
    fileID = fopen('landmarks_y.txt','w');
    fprintf(fileID, formatSpec, landmark_coords(:,2));
    fclose(fileID);
    
    LandmarksVisible = zeros(12,epochs/100);
    
else
    %BUILD IN LANDMARK SET UP FOR OTHER SIMULATION TYPES
    
    error('Landmark processing not programmed for these simulations');
    %     formatSpec = '%f\n';
    %
    %     fileID = fopen('landmarks_x.txt','w');
    %     fprintf(fileID, formatSpec, landmark_coords(:,1));
    %     fclose(fileID);
    %
    %     fileID = fopen('landmarks_y.txt','w');
    %     fprintf(fileID, formatSpec, landmark_coords(:,2));
    %     fclose(fileID);
end

disp('Landmarks set');


%% Allocating memory for normalisation operations
headings_vector = zeros((epochs/100)+1,1);
headings_matrix = zeros((epochs/100)+1,2);

HD_conjunctive_vector = zeros(1,conjunctive_cells);
HD_conjunctive_matrix = zeros(HD_cells,conjunctive_cells);

vis_conjunctive_vector = zeros(1,conjunctive_cells);
vis_conjunctive_matrix = zeros(num_landmarks,conjunctive_cells);

conjunctive_HD_vector = zeros(1,HD_cells);
conjunctive_HD_matrix = zeros(conjunctive_cells,HD_cells);

%WEIGHTS STATIC SO NOT NEEDED - LEFT THESE IN CASE LATER NEEDED
% IN CASE
% ADN_HD_vector = zeros(1,HD_cells);
% ADN_HD_matrix = zeros(ADN_cells,HD_cells);
%
% HD_ADN_vector = zeros(1,ADN_cells);
% HDA_ADN_matrix = zeros(HD_cells,ADN_cells);

%% Run a random walk path
%Set/Read Path Parameters
initial_position = [1,1];
initial_heading = [-0.5,0.5];
sigma_head_turns = 10;
step_scale_lower = 0.05;
step_scale_upper = 0.05;

if(strcmpi(model_type,'f')||strcmpi(model_type,'c'))
    [headings,positions] = one_compartment_walk_v2(0,epochs/100, initial_position,initial_heading,sigma_head_turns, step_scale_lower,step_scale_upper);
    disp('One compartment walk finished')
    
elseif(strcmpi(model_type,'b'))
    [headings,positions, odour_context] = two_compartment_walk(epochs/100, initial_position,initial_heading,sigma_head_turns, step_scale_lower,step_scale_upper);
    
    disp('Two compartment walk finished: ')
    vanilla_proportion = numel(odour_context(~odour_context))/(epochs/100);
    disp([num2str(vanilla_proportion*100),'% of time in vanilla']);
    lemon_proportion = numel(odour_context(odour_context>0.5))/(epochs/100);
    disp([num2str(lemon_proportion*100),'% of time in lemon']);
elseif(strcmpi(model_type,'t'))
    initial_position = [0.6,0.1]; %needs to be slightly different, as triangle not quite box
    [headings,positions, odour_context] = three_compartment_walk(epochs/100, initial_position,initial_heading,sigma_head_turns, step_scale_lower,step_scale_upper);
    disp('Three compartment walk finished')
    
else %this is the circular arena walk for visual scene simulation
    initial_position = [0,0]; %different as origin of circle
    sigma_head_turns = 15; %this chosen as works better for the circular arena simulations
    [headings,positions] = circular_arena_walk_v3(0,epochs/100, initial_position,initial_heading,sigma_head_turns, step_scale_lower,step_scale_upper);
    disp('Circular arena walk finished')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N.B. headings and positions give a value every %
% 100th timestep in the model! Addressed during  %
% the code, but be aware for future reference... %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get true HD (in deg) over course of simulation
[headings, ~, ~] = normalise_matrix(headings,'r',headings_vector,headings_matrix); %normalise each row of headings
trueHD = atan2(headings(:,2),headings(:,1));
trueHD(trueHD<0) = trueHD(trueHD<0) + (2*pi);
trueHD = trueHD * (180/pi); %convert to degrees

%% Allocate memory to store rates etc.
prev_rates_HD = zeros(HD_cells,1);
prev_activations_HD = zeros(HD_cells,1);

prev_rates_conjunctive = zeros(conjunctive_cells,1);
prev_activations_conjunctive = zeros(conjunctive_cells,1);

prev_rates_ADN = zeros(ADN_cells,1);
prev_activations_ADN = zeros(ADN_cells,1);

%Here's the rates over time - saving every 100th timestep to reduce memory
%usage

rates_HD_time = zeros(HD_cells,epochs/100);
rates_conjunctive_time = zeros(conjunctive_cells,epochs/100);
rates_vis_time = zeros(num_landmarks,epochs/100);
rates_ADN_time = zeros(ADN_cells,epochs/100);

disp('Network memory allocated')

%% Allocate memory to store weight change

%To keep track of weight change
HD_to_conjunctive_weight_change = zeros(epochs/100,1); %will be saved every 100 timesteps
vis_to_conjunctive_weight_change = zeros(epochs/100,1); %will be saved every 100 timesteps
conjunctive_to_HD_weight_change = zeros(epochs/100,1); %will be saved every 100 timesteps

%ADN TO RSC AND RSC TO ADN WEIGHTS STATIC - NO NEED TO RECORD WEIGHT
%CHANGES

%% Allocate memory for conduction buffers
BufferLength = conduction_delay/timestep_size;
circBuff_HD = zeros(HD_cells,BufferLength);
circBuff_conjunctive = zeros(conjunctive_cells,BufferLength);
circBuff_vis = zeros(num_landmarks,BufferLength);
circBuff_ADN = zeros(ADN_cells,BufferLength);

disp('Conduction delay buffers allocated');

%% Record the ADN increment each time HD changes
PI_increment_time = zeros(epochs/100,1);
true_increment_time = zeros(epochs/100,1);
ADN_PV_time = zeros(epochs/100,1);

%% Initialise weight profiles and connectivity
%set_favoured_view
[HD_PFDs] = set_PFDs(HD_cells);
CONJ_PFDs = [HD_PFDs,HD_PFDs]; %each half has own set of PFDs

%set HD to conjunctive weights
if(~strcmpi(model_type,'c'))
    [HD_to_conjunctive_weights,unconnected_conjunctives_HDpre] = initialise_weights('HD_to_conjunctive',HD_cells, conjunctive_cells, 'incomplete',50);
    fileID = fopen('HD_to_conjunctive_weights_initial.bdat','wb');
    fwrite(fileID,HD_to_conjunctive_weights,'float32'); %write to binary file for later use if needed
    initial_HD_to_conjunctive_weights = HD_to_conjunctive_weights; %record initial weights for plotting in script
    prev_HD_to_conjunctive_weights = HD_to_conjunctive_weights;
    
    fclose(fileID);
else
    [HD_to_conjunctive_weights, unconnected_conjunctives_HDpre] = pre_wire_weights('HD_to_conjunctive',HD_cells,conjunctive_cells,...
        HD_PFDs,CONJ_PFDs,HD_conjunctive_weight_sigma);
    fileID = fopen('HD_to_conjunctive_weights_initial.bdat','wb');
    fwrite(fileID,HD_to_conjunctive_weights,'float32'); %write to binary file for later use if needed
    initial_HD_to_conjunctive_weights = HD_to_conjunctive_weights; %record initial weights for plotting in script
    prev_HD_to_conjunctive_weights = HD_to_conjunctive_weights;
    fclose(fileID);
end

%set conjunctive to HD weights
if(~strcmpi(model_type,'c'))
    [conjunctive_to_HD_weights,unconnected_HD_conjunctivespre] = initialise_weights('conjunctive_to_HD',conjunctive_cells, HD_cells, 'complete',0);
    fileID = fopen('conjunctive_to_HD_weights_initial.bdat','wb');
    fwrite(fileID,conjunctive_to_HD_weights,'float32');%write to binary file for later use if needed
    initial_conjunctive_to_HD_weights = conjunctive_to_HD_weights; %record initial weights for plotting in script
    prev_conjunctive_to_HD_weights = conjunctive_to_HD_weights;
    fclose(fileID);
else
    [conjunctive_to_HD_weights, unconnected_HD_conjunctivespre] = pre_wire_weights('conjunctive_to_HD',conjunctive_cells,HD_cells,...
        CONJ_PFDs,HD_PFDs,conjunctive_HD_weight_sigma);
    fileID = fopen('conjunctive_to_HD_weights_initial.bdat','wb');
    fwrite(fileID,conjunctive_to_HD_weights,'float32');%write to binary file for later use if needed
    initial_conjunctive_to_HD_weights = conjunctive_to_HD_weights; %record initial weights for plotting in script
    prev_conjunctive_to_HD_weights = conjunctive_to_HD_weights;
    fclose(fileID);
end

%set vis to conjunctive weights
if(strcmpi(model_type,'f')||strcmpi(model_type,'v'))
    [vis_to_conjunctive_weights,unconnected_conjunctives_vispre] = initialise_weights('vis_to_conjunctive',num_landmarks, conjunctive_cells, 'complete',0);
    fileID = fopen('vis_to_conjunctive_weights_initial.bdat','wb');
    fwrite(fileID,vis_to_conjunctive_weights,'float32');
    %write to binary file for later use if needed
    initial_vis_to_conjunctive_weights = vis_to_conjunctive_weights; %record initial weights for plotting in script
    prev_vis_to_conjunctive_weights = vis_to_conjunctive_weights;
    fclose(fileID);
else
    [vis_to_conjunctive_weights, unconnected_conjunctives_vispre] = pre_wire_weights('vis_to_conjunctive',HD_cells,conjunctive_cells,...
        HD_PFDs,CONJ_PFDs,vis_conjunctive_weight_sigma);
    fileID = fopen('vis_to_conjunctive_weights_initial.bdat','wb');
    fwrite(fileID,vis_to_conjunctive_weights,'float32');
    %write to binary file for later use if needed
    initial_vis_to_conjunctive_weights = vis_to_conjunctive_weights; %record initial weights for plotting in script
    prev_vis_to_conjunctive_weights = vis_to_conjunctive_weights;
    fclose(fileID);
end

%set ADN to HD weights
[ADN_to_HD_weights,~] = pre_wire_weights('ADN_to_HD',ADN_cells,HD_cells,HD_PFDs,HD_PFDs,ADN_HD_weight_sigma);
fileID = fopen('ADN_to_HD_weights.bdat', 'wb');
fwrite(fileID,ADN_to_HD_weights,'float32');
fclose(fileID);

%set HD TO ADN weights
[HD_to_ADN_weights,~] = pre_wire_weights('HD_to_ADN',HD_cells,ADN_cells,HD_PFDs,HD_PFDs,HD_ADN_weight_sigma);
fileID = fopen('HD_to_ADN_weights.bdat', 'wb');
fwrite(fileID,HD_to_ADN_weights,'float32');
fclose(fileID);

%set ADN to ADN weights
[ADN_to_ADN_weights,~] = pre_wire_weights('ADN_to_ADN',ADN_cells,ADN_cells,HD_PFDs,HD_PFDs,ADN_ADN_weight_sigma);
fileID = fopen('ADN_to_ADN_weights.bdat', 'wb');
fwrite(fileID,ADN_to_ADN_weights,'float32');
fclose(fileID);

disp('Weights initialised and saved');

%% Converting percentages
PI_undershoot_percentage = PI_undershoot_percentage/100;
PI_noise_width = PI_noise_width/100;

disp('Press Any Key To Simulate');
pause();
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN SIMULATION
for timestep = 1:epochs
    
    disp(['Timestep: ',num2str(timestep)]);
    
    
    %% Set landmark cell firing
    if(strcmpi(model_type,'f')||strcmpi(model_type,'v')) %SO visual weights, set based on landmark visibility
        %landmark vectors
        landmark_vectors = landmark_coords - repmat(positions(ceil(timestep/100),:),size(landmark_coords,1),1);
        
        %angle between landmark vector and HD
        landmark_vectors = normalise_matrix(landmark_vectors,'r');
        
        curr_heading = headings(ceil(timestep/100),:);
        if strcmpi(model_type,'v')
            dot_prods = dot(landmark_vectors, curr_heading(ones(12,1),:),2);
        else
            dot_prods = dot(landmark_vectors, curr_heading(ones(num_landmarks,1),:),2);
        end
        
        %I think the landmark angles need to use the acos(dotproduct)
        %method here
        landmark_angles = acos(dot_prods); %angle between vectors a and b = acos(dot(a,b));
        
        if(strcmpi(model_type,'f'))
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
        
    else
        visdir = trueHD(ceil(timestep/100));
        
        if(strcmpi(model_type,'c'))
            visdir = visdir + conflict_size; %vision is offset for cue conflict
            visdir(visdir>360) = visdir(visdir>=360) - 360; %can't be more than 360
        else %bipolar/tripolar
            current_odour = odour_context(ceil(timestep/100));
            
            %bottom triangle = apple = 2
            %left triangle or compartment = vanilla = 1
            %right triangle or compartment = lemon = 0
            
            if current_odour>1.5     %apple
                %can only be in tripolar simulations
                visdir = visdir - 120;
            elseif current_odour >0.5 %lemon
                %+180 bipolar
                if(strcmpi(model_type,'b'))
                    visdir = visdir + 180;
                else %tripolar
                    visdir = visdir + 120;
                end
            end %vanilla is no change....
            
            visdir(visdir<0) = visdir(visdir<0) + 360; % make sure doesn't go past 360
            visdir(visdir>=360) = visdir(visdir>=360) - 360; %make sure doesn't go past 0
            
        end
        
        distance1 = abs((ones(1,HD_cells)*visdir) - HD_PFDs);
        distance2 = (ones(1,HD_cells)*360.0) - distance1;
        distance = min([distance1;distance2],[],1); %wraparound distance
        distance = (distance./vis_sigma).^2;
        
        %Guassian input, based on similarity of trueHD to PFD of each cell
        rates_vis = exp(-0.5.*distance);
        rates_vis = rates_vis'; %making correct dimension for later usage
    end
    
    
    %% Set path integration input to ADN
    if(timestep<2)
        PI_direction = trueHD(timestep); %first timestep ADN activity is initialised to true HD
    else
        %All other timesteps, ADN path integration input is noisy
        if trueHD(ceil(timestep/100))~=trueHD(ceil((timestep-1)/100)); %only a PI increment when trueHD changes
            
            %Initial increment is change in true HD
            PI_increment = atan2d(sind(trueHD(ceil(timestep/100))-trueHD(ceil((timestep-1)/100))),...
                cosd(trueHD(ceil(timestep/100))-trueHD(ceil((timestep-1)/100))));
            
            true_increment_time(ceil((timestep-1)/100)) = PI_increment; %record actual difference in HDs
            
            %Add noise, drawn from Guassian probility centered on increment
            %minus a percentage of increment, with sigma/width a certain
            %percentage of increment
            
            mean = PI_increment - (PI_undershoot_percentage*PI_increment);
            sigma = PI_noise_width*abs(PI_increment);
            
            PI_increment = normrnd(mean,sigma);
            
            PI_increment_time(ceil((timestep-1)/100)) = PI_increment;
            
            %Calculate ADN PVector from last timestep
            ADN_PV = calculate_PVector(prev_rates_ADN, HD_PFDs);
            ADN_PV_time(ceil((timestep-1)/100)) = ADN_PV;
            
            %New PI input direction is noisy path integration increment from
            %ADNPI
            
            PI_direction = PI_direction + PI_increment; %ADN_PV + PI_increment;
            
            if PI_direction > 360.0
                PI_direction = PI_direction - 360;
            end
            
            if PI_direction <= 0
                PI_direction = PI_direction + 360;
            end
            
        end
    end
    
    distance1 = abs((ones(1,HD_cells)*PI_direction) - HD_PFDs); %MUST BE ABS! %heading and position and true HD change every 100 timesteps....
    distance2 = (ones(1,HD_cells)*360.0) - distance1;
    distance = min([distance1;distance2],[],1); %wraparound distance
    distance = (distance./PI_sigma).^2;
    
    %Guassian input, based on similarity of trueHD to PFD of each cell
    PI_input = PI_strength .* exp(-0.5.*distance);
    PI_input = PI_input'; %making correct dimension for later usage
    
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Using Tony's trick to get 1000x500 prev_conjunctive rates from  %
    %   1000X1 prev_conjunctive_rates then multiply elementwise by    %
    %  weights to HD, then sum across HD dimension to give input for  %
    %  each HD cell - see normalise_matrix.m for details, or google   %
    %             tony's trick for more information                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %Old slower version
    %   conjunctive_buffer_rates = circBuff_conjunctive(:,end);
    %   expanded_conjunctive_buffer = conjunctive_buffer_rates(:,ones(HD_cells,1));
    %   conjunctive_input_matrix = expanded_conjunctive_buffer.*conjunctive_to_HD_weights;
    %   conjunctive_input = sum(conjunctive_input_matrix,1)';
    
    %Potentially faster bsxfun version
    conjunctive_input = sum(bsxfun(@times,conjunctive_to_HD_weights,circBuff_conjunctive(:,end)),1)';
    
    ADN_input = sum(bsxfun(@times,ADN_to_HD_weights,circBuff_ADN(:,end)),1)';
    
    activations_HD = (1.0 - HD_coefficient).*prev_activations_HD...
        + HD_coefficient.*(ADN_input.*ADN_HD_scale)...
        + HD_coefficient.*(conjunctive_input.*conjunctive_HD_scale)...
        - HD_coefficient.*(HD_inhibitory.*(1/HD_cells));
    
    
    %Conjunctive
    conjunctive_inhibitory = sum(prev_rates_conjunctive.*conjunctive_inhibition);
    
    %Old slower version
    %   HD_buffer_rates = circBuff_HD(:,end);
    %   expanded_HD_buffer = HD_buffer_rates(:,ones(conjunctive_cells,1));
    %   HD_input_matrix = expanded_HD_buffer.*HD_to_conjunctive_weights;
    %   HD_input = sum(HD_input_matrix,1)';
    
    %Potentially faster bsxfun version
    HD_input = sum(bsxfun(@times,HD_to_conjunctive_weights,circBuff_HD(:,end)),1)';
    
    
    
    %Old slower version
    %     vis_buffer_rates = circBuff_vis(:,end);
    %     expanded_vis_buffer = vis_buffer_rates(:,ones(conjunctive_cells,1));
    %     vis_input_matrix = expanded_vis_buffer.*vis_to_conjunctive_weights;
    %     vis_input = sum(vis_input_matrix,1)';
    
    
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
    
    
    if(timestep<=training_time/timestep_size) %don't update if past training phase
        %% Update weights
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    Weight update done by taking product of presynaptic    %
        %    and transposed postsynaptic, to yield outer product.   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %ORIGINAL SLOW VERSION
        %             HD_to_conjunctive_weights = prev_HD_to_conjunctive_weights + .... %Buffer Version
        %             (learning_rate*timestep_size).*(HD_buffer_rates*(prev_rates_conjunctive')); %prev_rates_conjunctive transposed to get outer product
        
        %FASTER BSXFUN VERSION
        if(LTD)
            conjunctive_avg_rate(1:conjunctive_cells) = nanmean(prev_rates_conjunctive);
            HD_to_conjunctive_weights = prev_HD_to_conjunctive_weights + ....
                (learning_rate*timestep_size).*bsxfun(@times,circBuff_HD(:,end),(prev_rates_conjunctive'-conjunctive_avg_rate));
        else
            HD_to_conjunctive_weights = prev_HD_to_conjunctive_weights + ....
                (learning_rate*timestep_size).*bsxfun(@times,circBuff_HD(:,end),(prev_rates_conjunctive'));
        end
        
        %Reset unconnected postsynaptic cell's receiving weights all to 0
        HD_to_conjunctive_weights(:,unconnected_conjunctives_HDpre) = 0;
        
        %VIS TO CONJUNCTIVE WEIGHTS WILL UPDATE FOR FULL SO, AND CUE
        %CONFLICT ONLY
        if(strcmpi(model_type,'f') || strcmpi(model_type,'c') || strcmpi(model_type,'v'))
            
            %ORIGINAL SLOW VERSION
            %             vis_to_conjunctive_weights = prev_vis_to_conjunctive_weights + .... %Buffer Version
            %                 (learning_rate*timestep_size).*(vis_buffer_rates*(prev_rates_conjunctive')); %prev_rates_conjunctive transposed to get outer product
            
            %FASTER BSXFUN VERSION
            if(LTD)
                vis_to_conjunctive_weights = prev_vis_to_conjunctive_weights + ....
                    (learning_rate*timestep_size).*bsxfun(@times,circBuff_vis(:,end),(prev_rates_conjunctive'-conjunctive_avg_rate));
            else
                vis_to_conjunctive_weights = prev_vis_to_conjunctive_weights + ....
                    (learning_rate*timestep_size).*bsxfun(@times,circBuff_vis(:,end),(prev_rates_conjunctive'));
            end
            
            %Reset unconnected postsynaptic cell's receiving weights all to 0
            vis_to_conjunctive_weights(:,unconnected_conjunctives_vispre) = 0;
        end
        
        %Conjunctive to HD
        
        %ORIGINAL SLOW VERSION
        %         conjunctive_to_HD_weights = prev_conjunctive_to_HD_weights + .... %Buffer Version
        %             (learning_rate*timestep_size).*(conjunctive_buffer_rates*(prev_rates_HD)'); %pre_rates_conjunctive transposed to get outer product
        
        %FASTER BSXFUN VERSION
        if(LTD)
            HD_avg_rate(1:HD_cells) = nanmean(prev_rates_HD);
            conjunctive_to_HD_weights = prev_conjunctive_to_HD_weights + ....
                (learning_rate*timestep_size).*bsxfun(@times,circBuff_conjunctive(:,end),(prev_rates_HD'-HD_avg_rate));
        else
            conjunctive_to_HD_weights = prev_conjunctive_to_HD_weights + ....
                (learning_rate*timestep_size).*bsxfun(@times,circBuff_conjunctive(:,end),(prev_rates_HD'));
        end
        
        conjunctive_to_HD_weights(:,unconnected_HD_conjunctivespre) = 0;
        
        %% Normalise weights
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Normalisation occurs by making all presynaptic weights for %
        %  a given postysnaptic cell have sum of squares equal to 1  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [HD_to_conjunctive_weights,HD_conjunctive_vector,HD_conjunctive_matrix] =...
            normalise_matrix(HD_to_conjunctive_weights,'c',HD_conjunctive_vector ,HD_conjunctive_matrix);
        
        [vis_to_conjunctive_weights,vis_conjunctive_vector,vis_conjunctive_matrix ] =...
            normalise_matrix(vis_to_conjunctive_weights,'c',vis_conjunctive_vector ,vis_conjunctive_matrix);
        
        [conjunctive_to_HD_weights,conjunctive_HD_vector,conjunctive_HD_matrix] =...
            normalise_matrix(conjunctive_to_HD_weights,'c',conjunctive_HD_vector ,conjunctive_HD_matrix);
        
        %Reset all unconnected postsynaptics receiving weights to 0, as a
        %failsafe after normalisation (technically not needed, but lol paranoia)
        
        HD_to_conjunctive_weights(:,unconnected_conjunctives_HDpre) = 0;
        vis_to_conjunctive_weights(:,unconnected_conjunctives_vispre) = 0;
        conjunctive_to_HD_weights(:,unconnected_HD_conjunctivespre) = 0;
        

    end %if training phase condition
    
    
    
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
        
        %Record weight changes
        
        %HD to Conjunctive
        weight_difference = abs(HD_to_conjunctive_weights - prev_HD_to_conjunctive_weights);
        HD_to_conjunctive_weight_change(timestep/100) = sum(weight_difference(:));
        
        %Vis to Conjunctive
        
        weight_difference = abs(vis_to_conjunctive_weights - prev_vis_to_conjunctive_weights);
        vis_to_conjunctive_weight_change(timestep/100) = sum(weight_difference(:));
        
        %Conjunctive to HD
        
        weight_difference = abs(conjunctive_to_HD_weights - prev_conjunctive_to_HD_weights);
        conjunctive_to_HD_weight_change(timestep/100) = sum(weight_difference(:));
        
        if(strcmpi(model_type,'v'))
            LandmarksVisible(:,timestep/100) = visible_landmarks;
        end
        
    end
    
    %Weights
    prev_HD_to_conjunctive_weights = HD_to_conjunctive_weights;
    prev_vis_to_conjunctive_weights = vis_to_conjunctive_weights;
    prev_conjunctive_to_HD_weights = conjunctive_to_HD_weights;
    
    
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
%% SAVING DATA
%buffers
fileID = fopen('VISBuffer.bdat','wb');
fwrite(fileID,circBuff_vis,'float32');
fclose(fileID);

fileID = fopen('HDBuffer.bdat','wb');
fwrite(fileID,circBuff_HD,'float32');
fclose(fileID);

fileID = fopen('CONJBuffer.bdat','wb');
fwrite(fileID,circBuff_conjunctive,'float32');
fclose(fileID);

fileID = fopen('ADNBuffer.bdat','wb');
fwrite(fileID,circBuff_ADN,'float32');
fclose(fileID);

disp('Conduction Buffers Saved');

%weight changes
fileID = fopen('VIStoCONJ_weight_change.bdat','wb');
fwrite(fileID,vis_to_conjunctive_weight_change,'float32');
fclose(fileID);

fileID = fopen('HDtoCONJ_weight_change.bdat','wb');
fwrite(fileID,HD_to_conjunctive_weight_change,'float32');
fclose(fileID);

fileID = fopen('CONJtoHD_weight_change.bdat','wb');
fwrite(fileID,conjunctive_to_HD_weight_change,'float32');
fclose(fileID);

%final weight profiles
fileID = fopen('HD_to_conjunctive_weights_final.bdat','wb');
fwrite(fileID,HD_to_conjunctive_weights,'float32');
fclose(fileID);

fileID = fopen('conjunctive_to_HD_weights_final.bdat','wb');
fwrite(fileID,conjunctive_to_HD_weights,'float32');
fclose(fileID);

fileID = fopen('vis_to_conjunctive_weights_final.bdat','wb');
fwrite(fileID,vis_to_conjunctive_weights,'float32');
fclose(fileID);

disp('Weight Data Saved');

%Rates
fileID = fopen('visRates.bdat','wb');
fwrite(fileID,rates_vis_time,'float32');
fclose(fileID);

fileID = fopen('HDRates.bdat','wb');
fwrite(fileID,rates_HD_time,'float32');
fclose(fileID);

fileID = fopen('conjunctiveRates.bdat','wb');
fwrite(fileID,rates_conjunctive_time,'float32');
fclose(fileID);

fileID = fopen('ADNRates.bdat','wb');
fwrite(fileID,rates_ADN_time,'float32');
fclose(fileID);

disp('Firing Rate Data Saved');

%save PVectors and PI increments
fileID = fopen('PI_increment.bdat','wb');
fwrite(fileID,PI_increment_time,'float32');
fclose(fileID);

fileID = fopen('true_increment.bdat','wb');
fwrite(fileID,true_increment_time,'float32');
fclose(fileID);

fileID = fopen('ADN_PV.bdat','wb');
fwrite(fileID,ADN_PV_time,'float32');
fclose(fileID);

%save landmark visibility if needed
if strcmpi(model_type,'v')
    fileID = fopen('LandmarkVisibility.bdat','wb');
    fwrite(fileID,LandmarksVisible,'float32');
    fclose(fileID);
end


figure()
subplot(2,1,1)
plot(true_increment_time,'r','Linewidth',2.0);
hold on
plot(PI_increment_time,'b','Linewidth',2.0);
title('Actual HD change vs. PI increment');
ylabel('Angle');
xlabel('Time');
subplot(2,1,2)
increment_difference = atan2d(sind(true_increment_time-PI_increment_time),cosd(true_increment_time-PI_increment_time));
increment_difference = abs(increment_difference);
increment_difference = increment_difference./abs(true_increment_time);
increment_difference = increment_difference.*100;
plot(increment_difference);
title('PI Error In Head Turns');
ylabel('Percentage');
xlabel('Time');
saveas(gcf,'_PI_inaccuracy', 'fig');
close(gcf);


% RSC_model_movie_v2(1,1,time,0.0001,500,250,250, 0, time,0);

%% PLOTTING DATA
%Weights

RSC_model_weights_v3('simulation',epochs,initial_vis_to_conjunctive_weights,vis_to_conjunctive_weights,initial_HD_to_conjunctive_weights...
    ,HD_to_conjunctive_weights,initial_conjunctive_to_HD_weights,conjunctive_to_HD_weights,vis_to_conjunctive_weight_change,...
    HD_to_conjunctive_weight_change,conjunctive_to_HD_weight_change, ADN_to_HD_weights, ADN_to_ADN_weights, HD_to_ADN_weights);


if(strcmpi(model_type,'f')||strcmpi(model_type,'c')||strcmpi(model_type,'v'))
    compartments = 1;
elseif(strcmpi(model_type,'b'))
    compartments = 2;
else
    compartments = 3;
end


%Rates
RSC_model_plot_v3('none','simulation',compartments,time, timestep_size, conjunctive_cells,num_landmarks, HD_cells, 'save',...
    rates_vis_time, rates_HD_time, rates_conjunctive_time,rates_ADN_time);

%% Disorientation test of ADN resetting via visual inputs
if(strcmpi(model_type,'v')) %only really worth doing for visual scene simulations
    disorientation_test(filename,positions(end-1,:),headings(end-1,:),1);
end


% %% Test of SO visual weights if needed
% if(strcmpi(model_type,'f'))
%     test_SO_RSC_v3(filename,positions(end-1,:),headings(end-1,:), 'absent',0);
% end
% if(strcmpi(model_type,'v'))
%     test_SO_RSC_v3(filename,positions(end-1,:),headings(end-1,:), 'absent',1);
% end


end
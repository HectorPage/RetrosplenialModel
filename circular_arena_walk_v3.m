function [headings,positions] = circular_arena_walk_v3(test_type,epochs,initial_position, first_heading,sigma,lower_scale_limit, upper_scale_limit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     GENERATING RANDOM WALK INSIDE CIRCULAR ARENA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Created by (and copyright held by Dr. Hector JI Page on 29/04/2016
%Modified by Dr. Hector JI Page on 23/01/2017
%Modified again by Dr. Hector JI Page on 03/02/2017

test_type = test_type(1); %taking first letter only to save on time with strcmpi

head_rotations = zeros(1,epochs+1);
head_rotations(1) = 0;

lower_scale_limit = abs(lower_scale_limit);   %set range of distance scaling factor - no <0 values
upper_scale_limit = abs(upper_scale_limit);

headings = zeros(epochs+1,2);
positions = zeros(epochs+1,2);

initial_heading = [0,0];
next_movement_direction = [0,0];

%Setting boundaries of envionment
%origin = [0,0];
circle_radius = 0.5; %radius in m

%%Reading in initial position

positions(1,:) = initial_position;  %recording the rat's start position

for counter = 1:epochs;
    disp(['Timestep: ',num2str(counter)]);
    
    %Randomly setting scaling factor for distance travelled.
    scaling = (upper_scale_limit-lower_scale_limit).*rand(1,1) + lower_scale_limit;
    
    %Setting initial position for this timestep
    initial_position = positions(counter,:);
    
    
    %%Setting and recording initial heading - could randomise this?
    if initial_heading(:) == 0;
        initial_heading = (first_heading/(norm(first_heading))) * scaling;    %[-0.5,0.5] originally
        headings(1,:) = initial_heading;
        movement_direction = initial_heading;
    else
        movement_direction = next_movement_direction;
    end
    
    movement_direction = movement_direction/(norm(movement_direction));
    
    %% Calculating final position and stopping rat from moving beyond maze edges
    
    movement_vector = movement_direction * scaling;
    final_position = initial_position + movement_vector;
    distance_from_origin = sqrt((0 - final_position(1)).^2   +   (0 - final_position(2)).^2); % distance calc.
    
    stopped_this_timestep = 0;
    while distance_from_origin>circle_radius %just keep reducing the movement
        movement_vector = movement_vector/2;
        final_position = initial_position + movement_vector;
        distance_from_origin = sqrt((0 - final_position(1)).^2   +   (0 - final_position(2)).^2); % distance calc.
        stopped_this_timestep = 1;
        
    end
    
    %% Recording final position
    
    positions(counter+1,:) = final_position;
    
    %% Calculating unrotated_final_heading
    
    movement_direction = movement_vector; %storing direction
    unrotated_final_heading = movement_direction/(norm(movement_direction));
    unrotated_final_heading = unrotated_final_heading/(norm(unrotated_final_heading));
    
    %% Calculating next motion direction (i.e. doing head rotation)
    
    %Here adding in a tendency to turn in the direction just turned
    if counter>1 %taking track of last turn taken
        if(position_rotation<180)
            right_turn = 1; %just turned R
        else
            right_turn = 0; %just turned L
        end
    end
    
    if stopped_this_timestep                    %Taking care of turning away from edges
        disp('Turning away from barrier');
        if(final_position(1)>0)
            if(final_position(2)>0)
                disp('North East');
            else
                disp('South East');
            end
        else
            if(final_position(2)>0)
                disp('North West');
            else
                disp('South West');
            end
        end
        LINE = [initial_position(1) initial_position(2)  unrotated_final_heading(1)  unrotated_final_heading(2)];
        CIRCLE = [0 0 0.5]; %central coordinates, and radius
        possible_points = intersectLineCircle(LINE, CIRCLE);
        
        %GET ANGLE BETWEEN MOVEMENT DIRECTION AND TANGENT
        %N.B. using the fact that radius and tangent are orthogonal to get
        %angle as 90 - angle between heading and radius
        
        %Angles below are correctd to account for quadrant (atan2d is
        %counter-clockwise by default)
        
        
        if(final_position(1)>0) %R hemisphere
            if(possible_points(1,1))>0 %choose intersection in R hemisphere
                intersection = possible_points(1,:);
            else
                intersection = possible_points(2,:);
            end
            
            %work out angle between tangent line and heading
            opp_mov_dir = -unrotated_final_heading;
            opp_mov_dir = opp_mov_dir/norm(opp_mov_dir);
            omd_azi = atan2(opp_mov_dir(2),opp_mov_dir(1));
            opp_rad = -intersection;
            opp_rad = opp_rad/norm(opp_rad);
            or_azi = atan2(opp_rad(2),opp_rad(1));
            
            angle = 90 - atan2d(sin(or_azi  - omd_azi),cos(or_azi  - omd_azi));
            
            if(final_position(2)>0) % TR quadrant
                %turn the correct direction
                
                %angle> 0 = turn L, otherwise turn R
                if(angle>0)
                    position_rotation = abs(normrnd(0,sigma)) * (pi/180); %turn left
                else
                    position_rotation = (0-abs(normrnd(0,sigma))) * (pi/180);  %turn right
                end
                
                
            else %otherwise BR quadrant
                %turn the correct direction
                
                %angle< 0 = turn L, otherwise turn R
                if(angle>0)
                    position_rotation = abs(normrnd(0,sigma)) * (pi/180); %turn left
                else
                    position_rotation = (0-abs(normrnd(0,sigma))) * (pi/180);  %turn right
                end
            end
            
        else %L hemisphere
            if(possible_points(1,1))<0 %choose intersection in R hemisphere
                intersection = possible_points(1,:);
            else
                intersection = possible_points(2,:);
            end
            
            %work out angle between tangent line and heading
            opp_mov_dir = -unrotated_final_heading;
            opp_mov_dir = opp_mov_dir/norm(opp_mov_dir);
            omd_azi = atan2(opp_mov_dir(2),opp_mov_dir(1));
            opp_rad = -intersection;
            opp_rad = opp_rad/norm(opp_rad);
            or_azi = atan2(opp_rad(2),opp_rad(1));
            
            angle = 90 - atan2d(sin(or_azi  - omd_azi),cos(or_azi  - omd_azi));
            if(final_position(2)>0) %TL quadrant
                
                %turn the correct direction
                %angle< 0 = turn L, otherwise turn R
                if(angle>90)
                    position_rotation = abs(normrnd(0,sigma)) * (pi/180); %turn left
                else
                    position_rotation = (0-abs(normrnd(0,sigma))) * (pi/180);  %turn right
                end
                
            else %otherwise BL quadrant
                %turn the correct direction
                
                %angle> 0 = turn L, otherwise turn R
                if(angle>90)
                    position_rotation = abs(normrnd(0,sigma)) * (pi/180); %turn left
                else
                    position_rotation = (0-abs(normrnd(0,sigma))) * (pi/180);  %turn right
                end
            end
        end
        
    else
        if counter>1
            continued_turn_probability = rand(1,1);
            if continued_turn_probability > 0.5 %50% chance of turning in the direction just taken
                if right_turn
                    position_rotation =  (0-abs(normrnd(0,sigma))) * (pi/180); %bias this towards the right
                else
                     position_rotation = abs(normrnd(0,sigma)) * (pi/180); %bias this towards the left
                end
            else
                position_rotation = normrnd(0,sigma) * (pi/180);  %normal turns
            end
        else
            position_rotation = normrnd(0,sigma) * (pi/180);  %normal turns  
        end
    end
    
    
    %% Correcting position rotation
    if position_rotation > (2*pi)
        position_rotation = (2*pi) - position_rotation;
    end%if position_rotation > 2pi
    
    head_rotations(counter+1) = position_rotation * (180/pi);
    
    %% Getting the new heading
    
    final_heading = [((cos(position_rotation)*unrotated_final_heading(1))-(sin(position_rotation)*unrotated_final_heading(2)))...
        , ((sin(position_rotation)*unrotated_final_heading(1))+(cos(position_rotation)*unrotated_final_heading(2)))];
    
    final_heading = final_heading/(norm(final_heading)) * scaling; %keep the headings the same length
    headings(counter+1,:) = final_heading;
    next_movement_direction = final_heading;
    
end %master loop

%% Working out head turns each timestep
trueHD = atan2(headings(:,2),headings(:,1));
trueHD(trueHD<0) = trueHD(trueHD<0) + (2*pi);
%trueHD = trueHD * (180/pi); %convert to degrees

PI_increment = zeros(numel(trueHD),1);

for idx = 2:numel(trueHD)
    PI_increment(idx) = atan2d(sin(trueHD(idx)-trueHD(idx-1)),...
        cos(trueHD(idx)-trueHD(idx-1)));
end

%% Creating a little movie to see how we're doing
% figure('doublebuffer','on','units','normalized','outerposition',[0 0 0.5 1]);
% for idx = 1:epochs
%     %plot a quiver
%     rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1]);
%     hold on
%     quiver(positions(idx,1),positions(idx,2),headings(idx,1),headings(idx,2),'b');
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%     ylim([-0.51 0.51]);
%     xlim([-0.51 0.51]);
%     axis equal
%     title(['Timestep: ',num2str(idx)]);
%     pause(0.01);
% end
%% Other summary plots characterising this walk

figure()
subplot(2,1,1)
histogram(PI_increment,-100:10:100)
xlim([-100,100])
xlabel('Rotation')
ylabel('Epochs')
set(gca,'xtick',-100:20:100)
set(gca,'Fontsize',24)
title('Head Turns Distribution')

subplot(2,1,2)
plot(PI_increment);
hold on
plot(head_rotations);
title('Head Turns Over Time')
savefig('_head_turn_stats');
close(gcf);

figure()
scatter(positions(1:end-1,1),positions(1:end-1,2),30,head_rotations(1:end-1),'filled');
colorbar;
caxis([-80 80]);
ylim([-0.51 0.51]);
xlim([-0.51 0.51]);
hold on;
rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Head Turn Locations');
savefig('_head_turn_locs');
close(gcf);


%% NOW PLOTTING HEADINGS
%THIS NEEDS TO BE EDITED
figure('Name','Headings','NumberTitle','off');
quiver(positions(:,1),positions(:,2),headings(:,1),headings(:,2),0.4,'b');
hold on;
plot(positions(1,1),positions(1,2),'X','MarkerSize',25,'MarkerFaceColor', 'r');
rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1]);
ylim([-0.51 0.51]);
xlim([-0.51 0.51]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
savefig('_path');
close(gcf);

%%
%Writing path data to text file for later

formatSpec = '%f\n';
if strcmpi(test_type,'s')
    fileID = fopen('positions_x_test.txt','w');
    fprintf(fileID, formatSpec, positions(:,1));
    fclose(fileID);
    
    fileID = fopen('positions_y_test.txt','w');
    fprintf(fileID, formatSpec, positions(:,2));
    fclose(fileID);
    
    fileID = fopen('headings_x_test.txt','w');
    fprintf(fileID, formatSpec, headings(:,1));
    fclose(fileID);
    
    fileID = fopen('headings_y_test.txt','w');
    fprintf(fileID, formatSpec, headings(:,2));
    fclose(fileID);
elseif strcmpi(test_type,'d')
    fileID = fopen('positions_x_disor.txt','w');
    fprintf(fileID, formatSpec, positions(:,1));
    fclose(fileID);
    
    fileID = fopen('positions_y_disor.txt','w');
    fprintf(fileID, formatSpec, positions(:,2));
    fclose(fileID);
    
    fileID = fopen('headings_x_disor.txt','w');
    fprintf(fileID, formatSpec, headings(:,1));
    fclose(fileID);
    
    fileID = fopen('headings_y_disor.txt','w');
    fprintf(fileID, formatSpec, headings(:,2));
    fclose(fileID);
else
    fileID = fopen('positions_x.txt','w');
    fprintf(fileID, formatSpec, positions(:,1));
    fclose(fileID);
    
    fileID = fopen('positions_y.txt','w');
    fprintf(fileID, formatSpec, positions(:,2));
    fclose(fileID);
    
    fileID = fopen('headings_x.txt','w');
    fprintf(fileID, formatSpec, headings(:,1));
    fclose(fileID);
    
    fileID = fopen('headings_y.txt','w');
    fprintf(fileID, formatSpec, headings(:,2));
    fclose(fileID);
end
end

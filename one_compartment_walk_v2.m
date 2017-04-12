function [headings,positions] = one_compartment_walk_v2(test_type,epochs,initial_position, first_heading,sigma,lower_scale_limit, upper_scale_limit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     GENERATING RANDOM WALK INSIDE SQUARE ARENA   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Created by (and copyright held by Dr. Hector JI Page on 29/04/2016
%Modified by Dr. Hector JI Page on 15/11/2016
%Modified again by Dr. Hector JI Page on 03/02/2017

%test types are 'SO_test' for a test without any PI input, or
%'disorientation' for a test with initial random ADN packet

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
NorthBound = 1.2; %ORIGINALLY 1.2
SouthBound = 0.0;
EastBound = 1.2;
WestBound = 0.0;

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
    
    %% Calculating final position
    
    final_position = initial_position + (movement_direction * scaling);
    
    %unrotated_final_heading = (final_position - initial_position); %final heading if rat had not rotated head
    
    pos_change = [initial_position(1) final_position(1); initial_position(2) final_position(2)]; %start and end point defining line of pos change
    
    
    %% Now stopping rat from moving beyond maze edges
    
    %define walls as lines based on bounds
    north_wall = [WestBound EastBound; NorthBound NorthBound]; %[x1,x2;y1,y2];
    south_wall = [WestBound EastBound; SouthBound SouthBound];
    east_wall = [EastBound EastBound; SouthBound NorthBound];
    west_wall = [WestBound WestBound; SouthBound NorthBound];
    
    %work out if pos_change crosses one of the walls
    North_intersection = InterX(pos_change,north_wall);
    South_intersection = InterX(pos_change,south_wall);
    East_intersection = InterX(pos_change,east_wall);
    West_intersection = InterX(pos_change,west_wall);
    
    %see if path has crossed any boundaries
    if(sum(North_intersection)>0)
        final_position = North_intersection' - ((North_intersection'-initial_position)/20); %stop on same direction but before wall
        
%        disp('Stopped At North Barrier');
    elseif(final_position(2)>NorthBound) %started on NorthBound, so no intersection
        final_position = initial_position;
    end
    
    if(sum(South_intersection)>0)
        final_position = South_intersection' - ((South_intersection'-initial_position)/20);
        
%        disp('Stopped At South Barrier');      
    elseif(final_position(2)<SouthBound) %started on SouthBound, so no intersection
        final_position = initial_position;
    end
    
    if(sum(East_intersection)>0)
        final_position = East_intersection' - ((East_intersection'-initial_position)/20);
        
%        disp('Stopped At East Barrier');
    elseif(final_position(1)>EastBound) %started on EastBound, so no intersection
        final_position = initial_position;
    end
    
    if(sum(West_intersection)>0)
        final_position = West_intersection' - ((West_intersection'-initial_position)/20);
        
%        disp('Stopped At West Barrier');
    elseif(final_position(1)<WestBound) %started on WestBound, so no intersection
        final_position = initial_position;
    end
    
    %%
    
    positions(counter+1,:) = final_position;
    
    %% Calculating unrotated_final_heading
    
    %below conditions allow for when the rat hasn't moved from same position
%     if ~any(final_position - initial_position)
        unrotated_final_heading = movement_direction/(norm(movement_direction));
%     else
%         unrotated_final_heading = (final_position - initial_position); %final heading if rat had not rotated head
%     end
    
    unrotated_final_heading = unrotated_final_heading/(norm(unrotated_final_heading));% * scaling; 
    
    
    %% Calculating next motion direction (i.e. head rotation)
    
    %%Taking care of turning away from edges
    
    if(initial_position(2) > NorthBound - scaling && unrotated_final_heading(2) > 0) %within a step of North, heading towards North
        %disp('Turned At North Barrier');
        if(unrotated_final_heading(1) > 0) %heading right
            position_rotation = (0-abs(normrnd(0,sigma))) * (pi/180);  %turn right
        else                               %heading left
            position_rotation = abs(normrnd(0,sigma)) * (pi/180); %turn left            
        end
    elseif(initial_position(2) < SouthBound+scaling && unrotated_final_heading(2) < 0) %within a step of South, heading towards South
        %disp('Turned At South Barrier');
        if(unrotated_final_heading(1) > 0) %heading right
            position_rotation = abs(normrnd(0,sigma)) * (pi/180); %turn left
        else                               %heading left
            position_rotation = (0-abs(normrnd(0,sigma))) * (pi/180); %turn right
        end
    elseif(initial_position(1) > EastBound-scaling && unrotated_final_heading(1) > 0) %within a step of East, heading towards East
        %disp('Turned At East Barrier');
        if(unrotated_final_heading(2) > 0) %heading up
            position_rotation = abs(normrnd(0,sigma)) * (pi/180); %turn left
        else                               %heading down
            position_rotation = (0-abs(normrnd(0,sigma))) * (pi/180); %turn right
        end
    elseif(initial_position(1) < WestBound+scaling && unrotated_final_heading(1) < 0) %within a step of West, heading towards West
       %disp('Turned At West Barrier');
        if(unrotated_final_heading(2) > 0) %heading up
              position_rotation = (0-abs(normrnd(0,sigma))) * (pi/180);  %turn right
        else                               %heading left
            position_rotation = abs(normrnd(0,sigma)) * (pi/180); %turn left   
        end
    else   
        position_rotation = normrnd(0,sigma) * (pi/180);  %normal turns
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

%% Creating headings along path

trueHD = atan2d(headings(:,2),headings(:,1));
trueHD(trueHD<0) = trueHD(trueHD<0) + (2*pi);
%trueHD = trueHD * (180/pi); %convert to degrees

PI_increment = zeros(numel(trueHD),1);
for idx = 2:numel(trueHD)
 PI_increment(idx) = atan2d(sind(trueHD(idx)-trueHD(idx-1)),...
                cosd(trueHD(idx)-trueHD(idx-1)));
end
% figure()
% plot(head_rotations)
% hold on
% plot(PI_increment)

figure()
subplot(2,1,1)
histogram(PI_increment,-200:10:200)
xlim([-200,200])
xlabel('Rotation')
ylabel('Epochs')
set(gca,'xtick',-180:40:180)
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
ylim([SouthBound NorthBound]);
xlim([WestBound EastBound]);

x = [WestBound, EastBound, EastBound, WestBound, WestBound];
y = [SouthBound, SouthBound, NorthBound, NorthBound, SouthBound];
hold on;
plot(x, y, 'k-', 'LineWidth', 1.5);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Head Turn Locations');
savefig('_head_turn_locs');
close(gcf);

%%
%%%NOW PLOTTING HEADINGS%%%
figure('Name','Headings','NumberTitle','off');
quiver(positions(:,1),positions(:,2),headings(:,1),headings(:,2),0.4,'b');
hold on;
plot(positions(1,1),positions(1,2),'X','MarkerSize',25,'MarkerFaceColor', 'r');
ylim([SouthBound-0.5 NorthBound+0.5]);
xlim([WestBound-0.5 EastBound+0.5]);

x = [WestBound, EastBound, EastBound, WestBound, WestBound];
y = [SouthBound, SouthBound, NorthBound, NorthBound, SouthBound];
plot(x, y, 'k-', 'LineWidth', 1.5);
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

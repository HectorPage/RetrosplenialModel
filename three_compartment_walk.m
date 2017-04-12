function [headings,positions, odour_context] = three_compartment_walk(epochs,initial_position, first_heading,sigma, lower_scale_limit, upper_scale_limit, SO_test)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GENERATING RANDOM WALK INSIDE TRIANGLUAR ARENA  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Created by (and copyright held by) Hector Page on 29/04/2016

head_rotations = zeros(1,epochs+1);
head_rotations(1) = 0;

lower_scale_limit = abs(lower_scale_limit);   %set range of distance scaling factor - no <0 values
upper_scale_limit = abs(upper_scale_limit);

headings = zeros(epochs+1,2);
positions = zeros(epochs+1,2);

initial_heading = [0,0];
next_movement_direction = [0,0];

%Setting boundaries of envionment
%triangle has length of sides of 1.2, therefore apex is at x=0.6
%and y=sqrt(3)/2 * 1.2, latter of which is not a whole number

%boundaries and barriers in format [x1, x2; y1, y2];

height = sqrt(3)/2 * 1.2;
centroid_x = 0.6;
centroid_y = height/3;

%trig to get the x and y change for the door in the left and right
%barriers, given door length of 0.05
x_change = 0.05*cosd(30);
y_change = 0.05 *(cosd(60));

%Boundaries
LeftBound = [0,centroid_x; 0,height];
BottomBound = [0,1.2; 0,0];
RightBound =[centroid_x,1.2; height,0];

%Barriers
TopBarrier = [centroid_x,centroid_x; height, centroid_y+0.05];
LeftBarrier = [0,centroid_x-x_change; 0,centroid_y-y_change];
RightBarrier = [1.2,centroid_x+x_change; 0,centroid_y-y_change];


%Apple compartment
Apple_compartment_x = [0 centroid_x 1.2];
Apple_compartment_y = [0 centroid_y 0];

%Vanilla compartment
Vanilla_compartment_x = [0 centroid_x centroid_x];
Vanilla_compartment_y = [0 height centroid_y];

%Lemon compartment
Lemon_compartment_x = [1.2 centroid_x centroid_x];
Lemon_compartment_y = [0 height centroid_y];

%Places to escape from corners to

vanilla_escape = [0.3, height/2];
lemon_escape = [0.9, height/2];
apple_escape = [0.6, centroid_y/2];


%%Reading in initial position

positions(1,:) = initial_position;  %recording the rat's start position

%initialising position rotation
position_rotation = 0;

%Odour contexts:
%bottom triangle = apple = 2
%left triangle = vanilla = 1
%right triangle = lemon = 0

odour_context = zeros(epochs,1);


for counter = 1:epochs;
    
    
%disp(['Path epoch: ',num2str(counter)]);
    %Randomly setting scaling factor for distance travelled.
    scaling = (upper_scale_limit-lower_scale_limit).*rand(1,1) + lower_scale_limit;
    
    %Setting initial position for this timestep
    initial_position = positions(counter,:);
    
    %Recording which odour context rat is currently in (assuming edge
    %never touched, and direct center never experienced
    [in_apple, onA] = inpolygon(positions(counter,1),positions(counter,2),Apple_compartment_x,Apple_compartment_y);
    [in_vanilla, onV] = inpolygon(positions(counter,1),positions(counter,2),Vanilla_compartment_x,Vanilla_compartment_y);
    [in_lemon,onL] =  inpolygon(positions(counter,1),positions(counter,2),Lemon_compartment_x,Lemon_compartment_y);
    
    if(in_apple)
        odour_context(counter) = 2;
    elseif(in_vanilla)
        odour_context(counter) = 1;
    elseif(in_lemon)
        odour_context(counter) = 0;
    else
        error('Compartment unknown!'); %if there's a failing of compartment detector
    end
    
    %% Setting and recording initial heading
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
    
    %% This section is to stop the rat hitting the barriers
    
    pos_change = [initial_position(1) final_position(1); initial_position(2) final_position(2)];
    
    TopBarr_intersection = InterX(pos_change,TopBarrier);
    LeftBarr_intersection = InterX(pos_change,LeftBarrier);
    RightBarr_intersection = InterX(pos_change,RightBarrier);
    
    %what about if rat hits top barrier coming from Apple?
    if(sum(TopBarr_intersection)>0) %If rat hits Top barrier
        if(odour_context(counter)==1) %if Vanilla
            final_position(1) = TopBarr_intersection(1) - ((scaling/300).*rand(1,1));
            final_position(2) = TopBarr_intersection(2);
        else %if Lemon
            final_position(1) = TopBarr_intersection(1) + ((scaling/300).*rand(1,1));
            final_position(2) = TopBarr_intersection(2);
            
        end
        %disp('Avoided Top Barrier');
    end
    
    %what about if rat hits left Barrier from Lemon?
    if(sum(LeftBarr_intersection)>0) %If rat hits Left barrier
        if(odour_context(counter)==1) %if Vanilla
            final_position(1) = LeftBarr_intersection(1);% - ((scaling/30).*rand(1,1));
            final_position(2) = LeftBarr_intersection(2) + ((scaling/300).*rand(1,1));
        elseif(odour_context(counter)==2) %if Apple
            final_position(1) = LeftBarr_intersection(1);% + ((scaling/30).*rand(1,1));
            final_position(2) = LeftBarr_intersection(2) - ((scaling/300).*rand(1,1));
            
        end
        
        %disp('Avoided Left Barrier');
    end
    
    %what about if rat hits Right Barrier from Vanilla?
    if(sum(RightBarr_intersection)>0) %If rat hits Right barrier
        if(odour_context(counter)<1) %if Lemon
            final_position(1) = RightBarr_intersection(1);% + ((scaling/30).*rand(1,1));
            final_position(2) = RightBarr_intersection(2) + ((scaling/300).*rand(1,1));
        elseif(odour_context(counter)==2) %if Apple
            final_position(1) = RightBarr_intersection(1);% - ((scaling/30).*rand(1,1));
            final_position(2) = RightBarr_intersection(2) - ((scaling/300).*rand(1,1));
            
        end
        
        %disp('Avoided Right Barrier');
    end
    
    %% Now stopping rat from moving beyond maze edges
    
    
    pos_change = [initial_position(1) final_position(1); initial_position(2) final_position(2)];
    
    left_intersection = InterX(pos_change,LeftBound);
    right_intersection = InterX(pos_change,RightBound);
    bottom_intersection = InterX(pos_change,BottomBound);
    
    
    if(sum(left_intersection)>0)         %LEFT BOUND HIT
        %only do this is exiting via vanilla, otherwise, an interior wall
        %is also blocking, so continue to use that
        if(in_vanilla)
            final_position(1) = left_intersection(1) + ((scaling/300).*rand(1,1));
            final_position(2) = left_intersection(2) - ((scaling/300).*rand(1,1));
        end
        
        %disp('Avoided Left Boundary');
              
        
    elseif(sum(right_intersection)>0)     %RIGHT BOUND HIT
        %only do this is exiting via lemon, otherwise, an interior wall
        %is also blocking, so continue to use that
        if(in_lemon)
            final_position(1) = right_intersection(1) - ((scaling/300).*rand(1,1));
            final_position(2) = right_intersection(2) - ((scaling/300).*rand(1,1));
        end
        %disp('Avoided Right Boundary');
    elseif(sum(bottom_intersection)>0)     %BOTTOM BOUND HIT
        %only do this if exiting via apple, otherwise, an interior wall
        %is also blocking, so continue to use that
        if(in_apple)
            final_position(1) = bottom_intersection(1);
            final_position(2) = bottom_intersection(2) + ((scaling/300).*rand(1,1));
        end
        %disp('Avoided Bottom Boundary');
    end%avoiding boundaries code
    
    
    %% NOW DOING ROTATIONS
    
    positions(counter+1,:) = final_position;
    
    %%Calculating unrotated_final_heading
    unrotated_final_heading = (final_position - initial_position); %final heading if rat had not rotated head
    unrotated_final_heading = unrotated_final_heading/(norm(unrotated_final_heading)) * scaling;
    
    %%Calculating next motion direction
    %To get the NEXT direction of movement:  do rotation by random amount
    
    %% Taking care of turning away from edges
    if(sum(left_intersection)>0)   %hit the Left Bound this timestep
        
        if(headings(counter,2)>=0 && headings(counter,1)<=-0.5)   %heading towards NNW
            position_rotation = normrnd(45,sigma) * (pi/180);
        elseif(headings(counter,2)>0 && headings(counter,1)<0 && headings(counter,1)>-0.5)%heading towards WNW
            position_rotation = normrnd(180,sigma) * (pi/180);
        elseif(headings(counter,2)<0 && headings(counter,1)<=-0.5) %heading WSW
            position_rotation = normrnd(-45,sigma) * (pi/180);
        end
        
    elseif(sum(right_intersection)>0)   %hit the Right Bound this timestep
        
        if(headings(counter,2)>=0 && headings(counter,1)<=0.5 && headings(counter,1)>0)   %heading towards NNE
            position_rotation = normrnd(-45,sigma) * (pi/180);
        elseif(headings(counter,2)>0 && headings(counter,1)>0 && headings(counter,1)>0.5)%heading towards ENE
            position_rotation = normrnd(180,sigma) * (pi/180);
        elseif(headings(counter,2)<0 && headings(counter,1)>=0.5) %heading ESE
            position_rotation = normrnd(45,sigma) * (pi/180);
        end
        
        
    elseif(sum(bottom_intersection)>0)   %hit the Bottom Bound this timestep
        
        if(headings(counter,1)<0)   %heading towards SW
            position_rotation = normrnd(60,sigma) * (pi/180);
        elseif(headings(counter,1)>0)%heading towards SE
            position_rotation = normrnd(-60,sigma) * (pi/180);
        end
        
    else
        %%Normal turning when not near an edge
        position_rotation = normrnd(0,sigma) * (pi/180);   %radian mode please
    end%turning from edges code
    
    
    
    %% This section is to rotate rat away from the barriers
    
    if(sum(TopBarr_intersection)>0) %If rat hits Top barrier - what about if from Apple?
        if(odour_context(counter)==1) %coming from Vanilla
            if(headings(counter,2)>0) %heading north
                position_rotation = normrnd(-45,sigma) * (pi/180);
            else%heading south
                position_rotation = normrnd(45,sigma) * (pi/180);
            end
        else %coming from Lemon
            if(headings(counter,2)>0) %heading north
                position_rotation = normrnd(45,sigma) * (pi/180);
            else%heading south
                position_rotation = normrnd(-45,sigma) * (pi/180);
            end
        end
    end %Top Barrier condition
    
    %what about if rat hits left Barrier from Lemon?
    if(sum(LeftBarr_intersection)>0) %If rat hits Left barrier
        if(odour_context(counter)==1) %if Vanilla, y component must be <0
            if(headings(counter,1)>=-0.5 && headings(counter,1)<0) %heading SSW
                position_rotation = normrnd(45,sigma) * (pi/180);
            elseif(headings(counter,1)<0.5) %heading SSE
                position_rotation = normrnd(180,sigma) * (pi/180);
            elseif(headings(counter,1)>=0.5) %heading ESE
                position_rotation = normrnd(-45,sigma) * (pi/180);
            end
        elseif(odour_context(counter)==2) %if Apple, y component must be >0
            if(headings(counter,1)<=-0.5)   %heading towards WNW
                position_rotation = normrnd(-45,sigma) * (pi/180);
            elseif(headings(counter,1)<0)%heading towards NNW
                position_rotation = normrnd(45,sigma) * (pi/180);
            end
        end
    end %Left Barrier condition
    
    %what about if rat hits Right Barrier from Vanilla?
    if(sum(RightBarr_intersection)>0) %If rat hits Right barrier
        if(odour_context(counter)<1) %if Lemon
            if(headings(counter,1)>=-0.5 && headings(counter,1)<0)%heading SSW
                position_rotation = normrnd(180,sigma) * (pi/180);
            elseif(headings(counter,1)<-0.5) %heading WSW
                position_rotation = normrnd(45,sigma) * (pi/180);
            elseif(headings(counter,1)>0 && headings(counter,1)<0.5) %heading SSE
                position_rotation = normrnd(-45,sigma) * (pi/180);
            end
        elseif(odour_context(counter)==2) %if Apple
            if(headings(counter,1)>=0.5)   %heading towards ENE
                position_rotation = normrnd(45,sigma) * (pi/180);
            elseif(headings(counter,1)>0)%heading towards NNE
                position_rotation = normrnd(-45,sigma) * (pi/180);
            end
        end
    end %Right Barrier condition
    
    
    if position_rotation > pi
        position_rotation = pi - position_rotation;
    end%if position_rotation > pi
    
    head_rotations(counter+1) = position_rotation;
    
    final_heading = [((cos(position_rotation)*unrotated_final_heading(1))-(sin(position_rotation)*unrotated_final_heading(2)))...
        , ((sin(position_rotation)*unrotated_final_heading(1))+(cos(position_rotation)*unrotated_final_heading(2)))];
    
    
     %% ADDING IN EXTRA CODE TO TEND TOWARDS THE DOOR
    
    door_probability = rand(1,1);
    
    if door_probability>0.3 % 70% probability of heading towards door with 50%
        if(final_position(2)>=centroid_y-0.15 && final_position(2)<=centroid_y+0.15) %if within vertical "gravity field" of door
            if(final_position(1)>=centroid_x-0.15 && final_position(1)<=centroid_x+0.15) %if within horizontal "gravity field" of door
                
                %find vector from final position to
                %centroid and work out angle between this and current
                %heading
                GV = [centroid_x-final_position(1),centroid_y-final_position(2)]; %goal vector from position
                heading_GV_angle = atan2d((final_heading(1)/norm(final_heading))*(GV(2)/norm(GV))-(final_heading(2)/norm(final_heading))*(GV(1)/norm(GV)),...
                    (final_heading(1)/norm(final_heading))*(GV(1)/norm(GV))+(final_heading(2)/norm(final_heading))*(GV(2)/norm(GV)));
                
                if(abs(heading_GV_angle)<=90) %if heading roughly towards the door - can change this threshold
                    door_position = [centroid_x+(-0.01+(0.01*rand(1,1))), centroid_y+(-0.01+(0.01*rand(1,1)))]; %the middle, jittered location
                    final_heading = door_position-final_position;
                end
                
            end
            
        end
    end
    
    
    %% ADDING IN EXTRA CODE TO STOP THE RAT GETTING TRAPPED IN CORNERS
    
    corner_escape_probability = rand(1,1);
    
    if corner_escape_probability >0.3 % 70% change of attempt to escape corner
        if(odour_context(counter)==1) %vanilla
            
            if(final_position(2)>=height-0.15)
                %find vector from final position to
                %escape and work out angle between this and current
                %heading
                EV = [vanilla_escape(1)-final_position(1),vanilla_escape(2)-final_position(2)]; %escape vector from position
                heading_EV_angle = atan2d((final_heading(1)/norm(final_heading))*(EV(2)/norm(EV))-(final_heading(2)/norm(final_heading))*(EV(1)/norm(EV)),...
                    (final_heading(1)/norm(final_heading))*(EV(1)/norm(EV))+(final_heading(2)/norm(final_heading))*(EV(2)/norm(EV)));
                
                if(abs(heading_EV_angle)<=120) %if heading roughly towards the escape - can vary angle
                    escape_position = [vanilla_escape(1)+(-0.001+(0.001*rand(1,1))), vanilla_escape(2)+(-0.001+(0.001*rand(1,1)))]; %jittered location
                    final_heading = escape_position-final_position;
                end
            end
            
            if(final_position(2)<=0.15 && final_position(1)<=0.15)
                %find vector from current position to
                %escape and work out angle between this and current
                %heading
                EV = [vanilla_escape(1)-final_position(1),vanilla_escape(2)-final_position(2)]; %escape vector from position
                heading_EV_angle = atan2d((final_heading(1)/norm(final_heading))*(EV(2)/norm(EV))-(final_heading(2)/norm(final_heading))*(EV(1)/norm(EV)),...
                    (final_heading(1)/norm(final_heading))*(EV(1)/norm(EV))+(final_heading(2)/norm(final_heading))*(EV(2)/norm(EV)));
                
                if(abs(heading_EV_angle)<=120) %if heading roughly towards the escape - can vary angle
                    escape_position = [vanilla_escape(1)+(-0.001+(0.001*rand(1,1))), vanilla_escape(2)+(-0.001+(0.001*rand(1,1)))]; %jittered location
                    final_heading = escape_position-final_position;
                end
            end
        end
        
        if(odour_context(counter)==0) %lemon
            if(final_position(2)>=height-0.15)
                %find vector from current position to
                %escape and work out angle between this and current
                %heading
                EV = [lemon_escape(1)-final_position(1),lemon_escape(2)-final_position(2)]; %escape vector from position
                heading_EV_angle = atan2d((final_heading(1)/norm(final_heading))*(EV(2)/norm(EV))-(final_heading(2)/norm(final_heading))*(EV(1)/norm(EV)),...
                    (final_heading(1)/norm(final_heading))*(EV(1)/norm(EV))+(final_heading(2)/norm(final_heading))*(EV(2)/norm(EV)));
                
                if(abs(heading_EV_angle)<=120) %if heading roughly towards the escape - can vary angle
                    escape_position = [lemon_escape(1)+(-0.001+(0.001*rand(1,1))), lemon_escape(2)+(-0.001+(0.001*rand(1,1)))]; %jittered location
                    final_heading = escape_position-final_position;
                end
            end
            
            if(final_position(2)<=0.15 && final_position(1)>=1.2-0.15)
                %find vector from current position to
                %escape and work out angle between this and current
                %heading
                EV = [lemon_escape(1)-final_position(1),lemon_escape(2)-final_position(2)]; %escape vector from position
                heading_EV_angle = atan2d((final_heading(1)/norm(final_heading))*(EV(2)/norm(EV))-(final_heading(2)/norm(final_heading))*(EV(1)/norm(EV)),...
                    (final_heading(1)/norm(final_heading))*(EV(1)/norm(EV))+(final_heading(2)/norm(final_heading))*(EV(2)/norm(EV)));
                
                if(abs(heading_EV_angle)<=120) %if heading roughly towards the escape - can vary angle
                    escape_position = [lemon_escape(1)+(-0.001+(0.001*rand(1,1))), lemon_escape(2)+(-0.001+(0.001*rand(1,1)))]; %jittered location
                    final_heading = escape_position-final_position;
                end
            end
        end
        
        
        if(odour_context(counter)==2) %apple
            
            if(final_position(2)<=0.15 && final_position(1)>=1.2-0.15)
                %find vector from current position to
                %escape and work out angle between this and current
                %heading
                EV = [apple_escape(1)-final_position(1),apple_escape(2)-final_position(2)]; %escape vector from position
                heading_EV_angle = atan2d((final_heading(1)/norm(final_heading))*(EV(2)/norm(EV))-(final_heading(2)/norm(final_heading))*(EV(1)/norm(EV)),...
                    (final_heading(1)/norm(final_heading))*(EV(1)/norm(EV))+(final_heading(2)/norm(final_heading))*(EV(2)/norm(EV)));
                
                if(abs(heading_EV_angle)<=120) %if heading roughly towards the escape - can vary angle
                    escape_position = [apple_escape(1)+(-0.001+(0.001*rand(1,1))), apple_escape(2)+(-0.001+(0.001*rand(1,1)))]; %jittered location
                    final_heading = escape_position-final_position;
                end
            end
            
            if(final_position(2)<=0.15 && final_position(1)<=0.15)
                %find vector from current position to
                %escape and work out angle between this and current
                %heading
                EV = [apple_escape(1)-final_position(1),apple_escape(2)-final_position(2)]; %escape vector from position
                heading_EV_angle = atan2d((final_heading(1)/norm(final_heading))*(EV(2)/norm(EV))-(final_heading(2)/norm(final_heading))*(EV(1)/norm(EV)),...
                    (final_heading(1)/norm(final_heading))*(EV(1)/norm(EV))+(final_heading(2)/norm(final_heading))*(EV(2)/norm(EV)));
                
                if(abs(heading_EV_angle)<=120) %if heading roughly towards the escape - can vary angle
                    escape_position = [apple_escape(1)+(-0.001+(0.001*rand(1,1))), apple_escape(2)+(-0.001+(0.001*rand(1,1)))]; %jittered location
                    final_heading = escape_position-final_position;
                end
            end
            
            
        end
        
        
    end %corner escape options
    
%      %% This section is to stop the rat from having performed an illegal movement based
%     
%     pos_change = [initial_position(1) final_position(1); initial_position(2) final_position(2)];
%     
%     left_intersection = InterX(pos_change,LeftBound);
%     right_intersection = InterX(pos_change,RightBound);
%     bottom_intersection = InterX(pos_change,BottomBound);
%     
%     TopBarr_intersection = InterX(pos_change,TopBarrier);
%     LeftBarr_intersection = InterX(pos_change,LeftBarrier);
%     RightBarr_intersection = InterX(pos_change,RightBarrier);
    
    %if any of these are true, make the final
    
    
    final_heading = final_heading/(norm(final_heading)) * scaling; %keep the headings the same length
    headings(counter+1,:) = final_heading;
    next_movement_direction = final_heading;
    
end %master loop


vanilla_proportion = numel(odour_context(~odour_context))/epochs;
disp([num2str(vanilla_proportion*100),'% of time in vanilla']);
lemon_proportion = numel(odour_context(odour_context>0.5 & odour_context<1.5))/epochs;
disp([num2str(lemon_proportion*100),'% of time in lemon']);
apple_proportion = numel(odour_context(odour_context>1.5))/epochs;
disp([num2str(apple_proportion*100),'% of time in apple']);



%Creating headings along path
for idx = 1:epochs-1
    headings(idx,:) = positions(idx+1,:) - positions(idx,:);
end

%%
%%%NOW PLOTTING HEADINGS%%%
figure('Name','Headings','NumberTitle','off');
quiver(positions(:,1),positions(:,2),headings(:,1),headings(:,2),0.4,'b');
hold on;
plot(positions(1,1),positions(1,2),'X','MarkerSize',25,'MarkerFaceColor', 'r');
ylim([0-0.5 height+0.5]);
xlim([0-0.5 1.2+0.5]);

x = [0,1.2, 0.6, 0];
y = [0, 0 , height, 0];
plot(x, y, 'k-', 'LineWidth', 1.5);

plot(TopBarrier(1,:), TopBarrier(2,:), 'k-', 'LineWidth', 1.5);
plot(LeftBarrier(1,:), LeftBarrier(2,:), 'k-', 'LineWidth', 1.5);
plot(RightBarrier(1,:), RightBarrier(2,:), 'k-', 'LineWidth', 1.5);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
savefig('_path');
close(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Writing path data to text file for later

formatSpec = '%f\n';

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

fileID = fopen('odour_context.txt','w');
fprintf(fileID, formatSpec, odour_context(:));
fclose(fileID);

end

function P = InterX(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve
%   together with any self-intersection points.
%
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')

%   Author : NS
%   Version: 3.0, 21 Sept. 2010

%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point.
%   Each factor of the 'C' arrays is essentially a matrix containing
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.

%...Argument checks and assignment of L2
error(nargchk(1,2,nargin));
if nargin == 1,
    L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
else
    L2 = varargin{1}; hF = @le;
end

%...Preliminary stuff
x1  = L1(1,:)';  x2 = L2(1,:);
y1  = L1(2,:)';  y2 = L2(2,:);
dx1 = diff(x1); dy1 = diff(y1);
dx2 = diff(x2); dy2 = diff(y2);

%...Determine 'signed distances'
S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);

C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

%...Obtain the segments where an intersection is expected
[i,j] = find(C1 & C2);
if isempty(i),P = zeros(2,0);return; end;

%...Transpose and prepare for output
i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0

%...Solve system of eqs to get the common points
P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
    dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';

    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end

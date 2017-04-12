function [headings,positions,odour_context] = two_compartment_walk(epochs,initial_position, first_heading,sigma, lower_scale_limit, upper_scale_limit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     GENERATING RANDOM WALK INSIDE SQUARE ARENA   %
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
NorthBound = 1.2;
SouthBound = 0.0;
EastBound = 1.2;
WestBound = 0.0;

%Adding in barrier elements
SouthBarrierStart = 0;
SouthBarrierEnd = 0.55;

NorthBarrierStart = 0.65;
NorthBarrierEnd = 1.2;

%%Reading in initial position

positions(1,:) = initial_position;  %recording the rat's start position

%initialising position rotation
position_rotation = 0;

odour_context = zeros(epochs,1); %This just says which half of the box rat is in, 1 for west, 0 for east.

for counter = 1:epochs;
    
    %Randomly setting scaling factor for distance travelled.
    scaling = (upper_scale_limit-lower_scale_limit).*rand(1,1) + lower_scale_limit;
    
    %Setting initial position for this timestep
    initial_position = positions(counter,:);
    
    %Recording whether in the west box
    if(initial_position(1) < 0.6)
        odour_context(counter) = 1;
    end
    
    
    %%Setting and recording initial heading - could randomise this?
    if initial_heading(:) == 0;
        initial_heading = (first_heading/(norm(first_heading))) * scaling;    %[-0.5,0.5] originally
        headings(1,:) = initial_heading;
        movement_direction = initial_heading;
    else
        movement_direction = next_movement_direction;
    end
    
    movement_direction = movement_direction/(norm(movement_direction));
    
    %%
    %%Calculating final position
    
    final_position = initial_position + (movement_direction * scaling);
    
    %unrotated_final_heading = (final_position - initial_position); %final heading if rat had not rotated head
    
    %%Now stopping rat from moving beyond maze edges%%
    %%
    if final_position(2) > NorthBound          %NORTH WALL
        if final_position(1) > EastBound       %NE corner
            final_position(1) = EastBound - ((scaling/2).*rand(1,1));
            final_position(2) = NorthBound - ((scaling/2).*rand(1,1)); %minus some limit maybe?
            
            %thestring = ['Stopped at NE corner: Timestep ',num2str(counter)];
            %disp(thestring);
        elseif final_position(1) < WestBound   %NW corner
            final_position(1) = WestBound + ((scaling/2).*rand(1,1));
            final_position(2) = NorthBound - ((scaling/2).*rand(1,1)); %minus some limit maybe?
            
            %thestring = ['Stopped at NW corner: Timestep ',num2str(counter)];
            %disp(thestring);
            
        else
            final_position(2) = NorthBound - ((scaling/2).*rand(1,1)); %minus some limit maybe?
            
            %thestring = ['Stopped at North wall: Timestep ',num2str(counter)];
            %disp(thestring);
        end
        
        
    elseif final_position(2) < SouthBound          %SOUTH WALL
        if final_position(1) > EastBound       %SE corner
            final_position(1) = EastBound - ((scaling/2).*rand(1,1));
            final_position(2) = SouthBound + ((scaling/2).*rand(1,1)); %minus some limit maybe?
            
            %thestring = ['Stopped at SE corner: Timestep ',num2str(counter)];
            %disp(thestring);
        elseif final_position(1) < WestBound   %NW corner
            final_position(1) = WestBound + ((scaling/2).*rand(1,1));
            final_position(2) = SouthBound + ((scaling/2).*rand(1,1)); %minus some limit maybe?
            
            %thestring = ['Stopped at SW corner: Timestep ',num2str(counter)];
            %disp(thestring);
            
        else
            final_position(2) = SouthBound + ((scaling/2).*rand(1,1)); %minus some limit maybe?
            
            %thestring = ['Stopped at South wall: Timestep ',num2str(counter)];
            %disp(thestring);
        end
        
    elseif final_position(1) > EastBound          %EAST WALL
        
        final_position(1) = EastBound - ((scaling/2).*rand(1,1)); %minus some limit maybe?
        
        %thestring = ['Stopped at East wall: Timestep ',num2str(counter)];
        %disp(thestring);
        
    elseif final_position(1) < WestBound          %WEST WALL
        
        final_position(1) = WestBound + ((scaling/2).*rand(1,1)); %minus some limit maybe?
        
        %thestring = ['Stopped at West wall: Timestep ',num2str(counter)];
        %disp(thestring);
    end%avoiding edges code
    
    %% THIS IS TO AVOID THE BARRIERS
    
    %% This section is to stop the rat hitting the barriers
    
    pos_change = [initial_position(1) final_position(1); initial_position(2) final_position(2)];
    south_barrier = [0.6 0.6; SouthBarrierStart SouthBarrierEnd];
    north_barrier = [0.6 0.6; NorthBarrierStart NorthBarrierEnd];
    
    south_intersection = InterX(pos_change,south_barrier);
    north_intersection = InterX(pos_change,north_barrier);
    
    
    if(sum(south_intersection)>0) %If rat hit South barrier
        if(initial_position(1)<0.6)
            final_position(1) = 0.59;
            final_position(2) = south_intersection(2);
        else
            final_position(1) = 0.61;
            final_position(2) = south_intersection(2);
            
        end
    end
    
    if(sum(north_intersection)>0) %If rat hit North barrier
        if(initial_position(1)<0.6)
            final_position(1) = 0.59;
            final_position(2) = north_intersection(2);
            
            
        else
            final_position(1) = 0.61;
            final_position(2) = north_intersection(2);
            
        end
    end
    
    
    %% NOW DOING ROTATIONS
    
    positions(counter+1,:) = final_position;
    
    %%Calculating unrotated_final_heading
    unrotated_final_heading = (final_position - initial_position); %final heading if rat had not rotated head
    unrotated_final_heading = unrotated_final_heading/(norm(unrotated_final_heading)) * scaling;
    
    %%Calculating next motion direction
    %To get the NEXT direction of movement:  do rotation by random amount
    
    %%Taking care of turning away from edges
    if(initial_position(2) > NorthBound-scaling) %within a step of North
        
        if(headings(counter,2)>0 && headings(counter,1)>0)   %heading towards northeast
            position_rotation = normrnd(-45,sigma) * (pi/180);
        elseif(headings(counter,2)>0 && headings(counter,1)<0)%heading towards northwest
            position_rotation = normrnd(45,sigma) * (pi/180);
        end
        
        if(initial_position(1) > EastBound-scaling)
            position_rotation = normrnd(180,sigma) * (pi/180);
        elseif(initial_position(1) < WestBound+scaling)
            position_rotation = normrnd(180, sigma) * (pi/180);
        end
        
    elseif(initial_position(2) < SouthBound+scaling) %within a step of South
        
        if(headings(counter,2)<0 && headings(counter,1)>0)   %heading towards southeast
            position_rotation = normrnd(45,sigma) * (pi/180);
        elseif(headings(counter,2)<0 && headings(counter,1)<0) %heading towards southwest
            position_rotation = normrnd(-45,sigma) * (pi/180);
        end
        
        if(initial_position(1) > EastBound-scaling)
            position_rotation = normrnd(180,sigma) * (pi/180);
        elseif(initial_position(1) < WestBound+scaling)
            position_rotation = normrnd(180, sigma) * (pi/180);
        end
        
        
    elseif(initial_position(1) > EastBound-scaling) %within a step of East
        if(headings(counter,1)>0 && headings(counter,2)>0)   %heading towards northeast
            position_rotation = normrnd(45,sigma) * (pi/180);
        elseif(headings(counter,1)>0 && headings(counter,2)<0) %heading towards southeast
            position_rotation = normrnd(-45,sigma) * (pi/180);
        end
        
        if(initial_position(2) > NorthBound-scaling)
            position_rotation = normrnd(180,sigma) * (pi/180);
        elseif(initial_position(2) < SouthBound+scaling)
            position_rotation = normrnd(180, sigma) * (pi/180);
        end
        
        
    elseif(initial_position(1) < WestBound+scaling) %within a step of West
        if(headings(counter,1)<0 && headings(counter,2)>0)   %heading towards northwest
            position_rotation = normrnd(-45,sigma) * (pi/180);
        elseif(headings(counter,1)<0 && headings(counter,2)<0) %heading towards southwest
            position_rotation = normrnd(45,sigma) * (pi/180);
        end
        
        if(initial_position(2) > NorthBound-scaling)
            position_rotation = normrnd(180,sigma) * (pi/180);
        elseif(initial_position(2) < SouthBound+scaling)
            position_rotation = normrnd(180, sigma) * (pi/180);
        end
    else
        %%Normal turning when not near an edge
        position_rotation = normrnd(0,sigma) * (pi/180);   %ditto, radian mode please
    end%turning from edges code
    
    %% This section is to rotate rat away from the barriers
    
    if(sum(south_intersection)>0) %If rat hit South barrier
        if(initial_position(1)<0.6)
            
            if(headings(counter,2) >0) %heading upwards
                position_rotation = normrnd(45,sigma) * (pi/180);
            else %heading downwards
                position_rotation = normrnd(-45,sigma) * (pi/180);
            end
        else
            
            if(headings(counter,2) >0) %heading upwards
                position_rotation = normrnd(-45,sigma) * (pi/180);
            else %heading downwards
                position_rotation = normrnd(45,sigma) * (pi/180);
            end
        end
    end
    
    if(sum(north_intersection)>0) %If rat hit North barrier
        if(initial_position(1)<0.6)
            
            if(headings(counter,2) >0) %heading upwards
                position_rotation = normrnd(-45,sigma) * (pi/180);
            else %heading downwards
                position_rotation = normrnd(45,sigma) * (pi/180);
            end
        else
            if(headings(counter,2) >0) %heading upwards
                position_rotation = normrnd(45,sigma) * (pi/180);
            else %heading downwards
                position_rotation = normrnd(-45,sigma) * (pi/180);
            end
        end
    end
    
    
    if position_rotation > pi
        position_rotation = pi - position_rotation;
    end%if position_rotation > pi
    
    head_rotations(counter+1) = position_rotation;
    
    final_heading = [((cos(position_rotation)*unrotated_final_heading(1))-(sin(position_rotation)*unrotated_final_heading(2)))...
        , ((sin(position_rotation)*unrotated_final_heading(1))+(cos(position_rotation)*unrotated_final_heading(2)))];
    
    
    
    %% ADDING IN EXTRA CODE TO TEND TOWARDS THE DOOR
    
    door_probability = rand(1,1);
    
    if door_probability>0.5 %probability of heading towards door with 50%
        if(initial_position(2)>=0.3 && initial_position(2)<=0.9) %if within vertical "gravity field" of door
            if(initial_position(1)>=0.4 && initial_position(1)<=0.6) %if within left "gravity field" of door
                
                if(unrotated_final_heading(1)>0 && abs(unrotated_final_heading(2))<0.5) %if heading roughly towards the door
                    
                    door_position = [0.6+(-0.01+(0.01*rand(1,1))), 0.6+(-0.01+(0.01*rand(1,1)))]; %the middle, jittered location
                    final_heading = door_position-initial_position;
                    
                end
                
            end
            
            if(initial_position(1)>0.6 && initial_position(1)<=0.8) %if within right "gravity field" of door
                
                if(unrotated_final_heading(1)<0 && abs(unrotated_final_heading(2))<0.5) %if heading roughly towards the door
                    
                    door_position = [0.6+(-0.01+(0.01*rand(1,1))), 0.6+(-0.01+(0.01*rand(1,1)))]; %the middle, jittered location
                    final_heading = door_position-initial_position;
                    
                end
            end
        end
    end
    
    
    final_heading = final_heading/(norm(final_heading)) * scaling; %keep the headings the same length
    headings(counter+1,:) = final_heading;
    next_movement_direction = final_heading;
    
end %master loop

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
ylim([SouthBound-0.5 NorthBound+0.5]);
xlim([WestBound-0.5 EastBound+0.5]);

x = [WestBound, EastBound, EastBound, WestBound, WestBound];
y = [SouthBound, SouthBound, NorthBound, NorthBound, SouthBound];
plot(x, y, 'k-', 'LineWidth', 1.5);

x1 = [0.6 0.6];
y1 = [0 0.55];
y2 = [0.65 1.2];

plot(x1, y1, 'k-', 'LineWidth', 1.5);
plot(x1, y2, 'k-', 'LineWidth', 1.5);
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

function RSC_model_movie_v3(compartments,fully_SO,seconds,tsize,conjunctive_cells,hd_cells,adn_cells,start_time, end_time,test_type,visual_scene)
%% MOVIE OF NETWORK BEHAVIOUR
test_type = test_type(1);
disp('Loading Data...')
%% SETTING UP TIMING INFO
steps = uint32(((seconds/tsize)/100)); %Dividing by 100, as new program saves every 100th timestep, to reduce filesize.

if ~start_time %can just set at 0 to go from start
    start_steps = 1;
else
    start_steps = uint32(((start_time/tsize)/100));
end
end_steps = uint32(((end_time/tsize)/100));

if ~fully_SO
    vis_cells = hd_cells; %this is typically true of my simulations, but be aware of this
else
    vis_cells = 2; %originally 12 of them
end

if visual_scene
    vis_cells = 12*6;
end


%% Getting True HD over course of simulation
if strcmpi(test_type,'s')
    fid = fopen('headings_x_test.txt', 'r');
    head_x = fscanf(fid ,'%f\n');
    fclose(fid);
    
    fid = fopen('headings_y_test.txt', 'r');
    head_y = fscanf(fid,'%f\n');
    fclose(fid);
elseif (strcmpi(test_type,'d'))
    fid = fopen('headings_x_disor.txt', 'r');
    head_x = fscanf(fid ,'%f\n');
    fclose(fid);
    
    fid = fopen('headings_y_disor.txt', 'r');
    head_y = fscanf(fid,'%f\n');
    fclose(fid);
else
    fid = fopen('headings_x.txt', 'r');
    head_x = fscanf(fid ,'%f\n');
    fclose(fid);
    
    fid = fopen('headings_y.txt', 'r');
    head_y = fscanf(fid,'%f\n');
    fclose(fid);
end

%Normalise all headings
head_lengths = sqrt((head_x.^2)+(head_y.^2));

head_x = head_x./head_lengths;
head_y = head_y./head_lengths;

%Work out true HD over course of simulation
trueHD = atan2(head_y,head_x);
trueHD(trueHD<0) = trueHD(trueHD<0) + (2*pi);
trueHD = trueHD * (180/pi);

%% Getting positions
if strcmpi(test_type,'s')
    fid = fopen('positions_x_test.txt', 'r');
    pos_x = fscanf(fid ,'%f\n');
    fclose(fid);
    
    fid = fopen('positions_y_test.txt', 'r');
    pos_y = fscanf(fid,'%f\n');
    fclose(fid);
elseif strcmpi(test_type,'d')
     fid = fopen('positions_x_disor.txt', 'r');
    pos_x = fscanf(fid ,'%f\n');
    fclose(fid);
    
    fid = fopen('positions_y_disor.txt', 'r');
    pos_y = fscanf(fid,'%f\n');
    fclose(fid);
else
    fid = fopen('positions_x.txt', 'r');
    pos_x = fscanf(fid ,'%f\n');
    fclose(fid);
    
    fid = fopen('positions_y.txt', 'r');
    pos_y = fscanf(fid,'%f\n');
    fclose(fid);
end

clear head_lengths


% %% Getting compartment
% if compartments>1.5
%     fid = fopen('odour_context.txt', 'r');
%     odour_context = fscanf(fid ,'%f\n');
%     fclose(fid);
% end

%% Getting landmark visibility
if visual_scene
    LandmarkVisibility = file_load(steps,12,'LandmarkVisibility.bdat'); %these are which landmarks are visible (different to vis cell firing for these simulations)
    if strcmpi(test_type,'s')
        LandmarkVisibility = file_load(steps,12,'LandmarkVisibility_test.bdat'); %these are which landmarks are visible (different to vis cell firing for these simulations)
    elseif(strcmpi(test_type,'d'))
         LandmarkVisibility = file_load(steps,12,'LandmarkVisibility_disor.bdat'); %these are which landmarks are visible (different to vis cell firing for these simulations)
    end
end

%% Getting rates of the network elements
if strcmpi(test_type,'s')
    ConjunctiveRates = file_load(steps, conjunctive_cells, 'conjunctiveRates_test.bdat')'; %these are transposed, as data output different dimensions to C codes
    HDRates = file_load(steps,hd_cells,'HDRates_test.bdat')';
    VisRates = file_load( steps,vis_cells, 'visRates_test.bdat')';
    ADNRates = file_load(steps,adn_cells,'ADNRates_test.bdat')';
elseif strcmpi(test_type,'d')
    ConjunctiveRates = file_load(steps, conjunctive_cells, 'conjunctiveRates_disor.bdat')'; %these are transposed, as data output different dimensions to C codes
    HDRates = file_load(steps,hd_cells,'HDRates_disor.bdat')';
    VisRates = file_load( steps,vis_cells, 'visRates_disor.bdat')';
    ADNRates = file_load(steps,adn_cells,'ADNRates_disor.bdat')';
else
    ConjunctiveRates = file_load(steps, conjunctive_cells, 'conjunctiveRates.bdat')'; %these are transposed, as data output different dimensions to C codes
    HDRates = file_load(steps,hd_cells,'HDRates.bdat')';
    VisRates = file_load( steps,vis_cells, 'visRates.bdat')';
    ADNRates = file_load(steps,adn_cells,'ADNRates.bdat')';
end

%% Getting directions (PVectors) of each network region

disp('Generating LSE and RSC Population Vectors....')
%Caluclating PVs

%favoured view needed
increment = 360/hd_cells;
favoured_view = (0:hd_cells-1)*increment;
favoured_view = favoured_view';

fav_view_sin = sind(favoured_view); %sine
fav_view_cos = cosd(favoured_view); %cosine

if ~fully_SO
    vector_1_vis = sum(bsxfun(@times,VisRates,fav_view_sin),1);
    vector_2_vis = sum(bsxfun(@times,VisRates,fav_view_cos),1);
    
    %Get atan2 (four-quadrant corrected) angles
    VisPV = atan2d(vector_1_vis,vector_2_vis);
    
    %Get in the right range (0-360)
    VisPV(VisPV<0) = VisPV(VisPV<0) + 360;
end

vector_1_rsc = sum(bsxfun(@times,HDRates,fav_view_sin),1);
vector_2_rsc = sum(bsxfun(@times,HDRates,fav_view_cos),1);

%Get atan2 (four-quadrant corrected) angles
RSCPV = atan2d(vector_1_rsc,vector_2_rsc);

%Get in the right range (0-360)
RSCPV(RSCPV<0) = RSCPV(RSCPV<0)  + 360;

vector_1_adn = sum(bsxfun(@times,ADNRates,fav_view_sin),1);
vector_2_adn = sum(bsxfun(@times,ADNRates,fav_view_cos),1);

%Get atan2 (four-quadrant corrected) angles
ADNPV = atan2d(vector_1_adn,vector_2_adn);

%Get in the right range (0-360)
ADNPV(ADNPV<0) = ADNPV(ADNPV<0)  + 360;


%% Setting up subplots for movie

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subplots panels for each type

vis_panels = [1 2 3 7 8 9];
rsc_panels = [4 5 6 10 11 12];
path_panels = [14 15 20 21];
ADN_panels = [16 17 18 22 23 24];
WC_panels = [25 26 27 31 32 33];
BC_panels = [28 29 30 34 35 36];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Visual Cells
figure('doublebuffer','on','units','normalized','outerposition',[0 0 0.5 1]);  %half-screen buffered figure
hold on;
subplot(6,6,vis_panels);
set(gca,'box','on');
axis equal
set(gca,'Xlim',[0 500]);
set(gca,'ylim',[0 1.1]);
axis manual;
grid off;
set(gca,'YTick',0:0.25:1);
set(gca,'XTick',0:125:500);
set(gca,'XTickLabel',0:125*increment:360);
title('VIS');

if fully_SO
    %Read in landmarks coordinates for N landmarks
    fid = fopen('landmarks_x.txt', 'r');
    landmark_x = fscanf(fid ,'%f\n');
    fclose(fid);
    
    fid = fopen('landmarks_y.txt', 'r');
    landmark_y = fscanf(fid,'%f\n');
    fclose(fid);
    
    landmark_size = 0.025; %for plotting purposes
    
    if visual_scene
        
        landmark_x(landmark_x>0) = landmark_x(landmark_x>0)-landmark_size;
        landmark_y(landmark_y>0) = landmark_y(landmark_y>0)-landmark_size;
        
    else
        
        landmark_x(landmark_x>1.1) = landmark_x(landmark_x>1.1)-landmark_size;
        landmark_y(landmark_y>1.1) = landmark_y(landmark_y>1.1)-landmark_size;
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RSCHD Cells
subplot(6,6,rsc_panels);
set(gca,'box','on');
axis equal
set(gca,'Xlim',[0 hd_cells]);
set(gca,'ylim',[0 1.1]);
axis manual;
grid off;
set(gca,'YTick',0:0.25:1);
set(gca,'XTick',0:hd_cells/4:hd_cells);
set(gca,'XTickLabel',0:(hd_cells/4)*increment:360);
title('RSC HD');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WC Cells
subplot(6,6,WC_panels);
set(gca,'box','on');
axis equal
set(gca,'Xlim',[0 conjunctive_cells]);
set(gca,'ylim',[0 1.1]);
axis manual;
grid off;
set(gca,'YTick',0:0.25:1);
set(gca,'XTick',0:conjunctive_cells/4:conjunctive_cells);
set(gca,'XTickLabel',0:(conjunctive_cells/4)*increment:360);
title('WC Flip');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BC Cells
subplot(6,6,BC_panels);
set(gca,'box','on');
axis equal
set(gca,'Xlim',[0 conjunctive_cells]);
set(gca,'ylim',[0 1.1]);
axis manual;
grid off;
set(gca,'YTick',0:0.25:1);
set(gca,'XTick',0:conjunctive_cells/4:conjunctive_cells);
set(gca,'XTickLabel',0:(conjunctive_cells/4)*increment:360);
title('BC Flip');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HD Plot
subplot(6,6,ADN_panels);
set(gca,'box','on');
axis equal
set(gca,'Xlim',[0 500]);
set(gca,'ylim',[0 1.1]);
axis manual;
grid off;
set(gca,'YTick',0:0.25:1);
set(gca,'XTick',0:hd_cells/4:hd_cells);
set(gca,'XTickLabel',0:(hd_cells/4)*increment:360);
title('ADN');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setting up info for drawing path

if compartments>2.5 %FOR TRIPOLAR SIMULATIONS
    height = sqrt(3)/2 * 1.2;
    centroid_x = 0.6;
    centroid_y = height/3;
    
    %trig to get the x and y change for the door in the left and right
    %barriers, given door length of 0.05
    x_change = 0.05*cosd(30);
    y_change = 0.05 *(cosd(60));
    
    %Barriers
    TopBarrier = [centroid_x,centroid_x; height, centroid_y+0.05];
    LeftBarrier = [0,centroid_x-x_change; 0,centroid_y-y_change];
    RightBarrier = [1.2,centroid_x+x_change; 0,centroid_y-y_change];
    
    
end

% For the rat triangle
headWidth = 5;
headLength = 5;
lineLength = 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Initial Frame
mov=getframe(1);
%% Central Plotting Script

for timestep = start_steps:end_steps
    
    % VIS
    subplot(6,6,vis_panels);
    if ~fully_SO %only plot vis activation if not fully SO, otherwise just put landmarks on path plot
        plot(VisRates(:,timestep),'b','Linewidth',2.0);
        hold on
        plot([VisPV(timestep)/increment VisPV(timestep)/increment],[0 1],'Color',[0 0 1],'Linewidth',2.0, 'Linestyle','--');
        plot([trueHD(timestep)/increment trueHD(timestep)/increment],[0 1],'Color',[0 0 0 0.5],'Linewidth',2.0);
        set(gca,'XTick',0:vis_cells/4:vis_cells);
        set(gca,'XTickLabel',0:(vis_cells/4)*increment:360);
        hold off
        title('VIS');
    end
    
    % RSC HD
    subplot(6,6,rsc_panels);
    plot(HDRates(:,timestep),'r','Linewidth',2.0);
    hold on
    if ~fully_SO
        plot([VisPV(timestep)/increment VisPV(timestep)/increment],[0 1],'Color',[0 0 1],'Linewidth',2.0, 'Linestyle','--');
    end
    plot([RSCPV(timestep)/increment RSCPV(timestep)/increment],[0 1],'Color',[1 0 0],'Linewidth',2.0, 'Linestyle',':');
    plot([trueHD(timestep)/increment trueHD(timestep)/increment],[0 1],'Color',[0 0 0 0.5],'Linewidth',2.0);
    plot([ADNPV(timestep)/increment ADNPV(timestep)/increment],[0 1],'Color','m','Linewidth',2.0, 'Linestyle','-');
    set(gca,'XTick',0:hd_cells/4:hd_cells);
    set(gca,'XTickLabel',0:(hd_cells/4)*increment:360);
    hold off
    title('RSC HD');
    
    if fully_SO
        % WC and BC distinction meaningless - plot whole conjunctive layer
        all_conjunctive = [WC_panels,BC_panels]; %plot entire layer
        conjunctive_increment = 360/conjunctive_cells;
        subplot(6,6,all_conjunctive);
        plot(ConjunctiveRates(:,timestep),'Color','k','Linewidth',2.0);
        set(gca,'XTick',0:(conjunctive_cells/4):conjunctive_cells);
        set(gca,'XTickLabel',0:(conjunctive_cells/4)*conjunctive_increment:360);
        xlim([0 conjunctive_cells]);
        ylim([0 1.1]);
        hold off
        title('Conjunctive Cells');
    else
        
        % WC Flips
        subplot(6,6,WC_panels);
        plot(ConjunctiveRates(1:vis_cells,timestep),'Color','k','Linewidth',2.0);
        hold on
        plot([VisPV(timestep)/increment VisPV(timestep)/increment],[0 1],'Color',[0 0 1],'Linewidth',2.0, 'Linestyle','--');
        plot([RSCPV(timestep)/increment RSCPV(timestep)/increment],[0 1],'Color',[1 0 0],'Linewidth',2.0, 'Linestyle',':');
        plot([trueHD(timestep)/increment trueHD(timestep)/increment],[0 1],'Color',[0 0 0 0.5],'Linewidth',2.0);
        set(gca,'XTick',0:hd_cells/4:hd_cells);
        set(gca,'XTickLabel',0:(hd_cells/4)*increment:360);
        hold off
        title('WC');
        
        % BC Flips
        subplot(6,6,BC_panels);
        plot(ConjunctiveRates(vis_cells+1:vis_cells*2,timestep),'Color','k','Linewidth',2.0);
        hold on
        plot([VisPV(timestep)/increment VisPV(timestep)/increment],[0 1],'Color',[0 0 1],'Linewidth',2.0, 'Linestyle','--');
        plot([RSCPV(timestep)/increment RSCPV(timestep)/increment],[0 1],'Color',[1 0 0],'Linewidth',2.0, 'Linestyle',':');
        plot([trueHD(timestep)/increment trueHD(timestep)/increment],[0 1],'Color',[0 0 0 0.5],'Linewidth',2.0);
        set(gca,'XTick',0:hd_cells/4:hd_cells);
        set(gca,'XTickLabel',0:(hd_cells/4)*increment:360);
        hold off
        title('BC');
    end
    
    % True HD
    subplot(6,6,ADN_panels);
    plot(ADNRates(:,timestep),'m','Linewidth',2.0);
    hold on
    if ~fully_SO
        plot([VisPV(timestep)/increment VisPV(timestep)/increment],[0 1],'Color',[0 0 1],'Linewidth',2.0, 'Linestyle','--');
    end
    
    plot([RSCPV(timestep)/increment RSCPV(timestep)/increment],[0 1],'Color',[1 0 0],'Linewidth',2.0, 'Linestyle',':');
    plot([ADNPV(timestep)/increment ADNPV(timestep)/increment],[0 1],'Color','m','Linewidth',2.0, 'Linestyle','-');
    plot([trueHD(timestep)/increment trueHD(timestep)/increment],[0 1],'Color',[0 0 0 0.5],'Linewidth',2.0);
    set(gca,'XTick',0:hd_cells/4:hd_cells);
    set(gca,'XTickLabel',0:(hd_cells/4)*increment:360);
    hold off
    title('ADN');
    
    % Path
    if fully_SO
        landmark_panels = [vis_panels,path_panels];
        subplot(6,6,landmark_panels);
    else
        subplot(6,6,path_panels);
    end
    
    if fully_SO
        if visual_scene
            hold on
            rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1]);
            
            for i=1:12 %looping over all the visible landmarks
                
                rectangle('position', [landmark_x(i) landmark_y(i) landmark_size landmark_size],...
                    'facecolor', ones(1,3)-LandmarkVisibility(timestep,i), 'edgecolor', 'k');
                
            end
            
            ylim([-0.51 0.51]);
            xlim([-0.51 0.51]);
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            title('LANDMARKS')
            
        else
            %Fully SO - show the box and the landmarks (black if visible, white
            %if not visible)
            hold on
            x = [0, 1.2, 1.2, 0];
            y = [0, 0 , 1.2, 1.2];
            plot(x, y, 'k-', 'LineWidth', 1.0);
            plot([0 0],[0 1.2], 'k-','LineWidth',1.0);
            for i=1:vis_cells,
                
                rectangle('position', [landmark_x(i) landmark_y(i) landmark_size landmark_size],...
                    'facecolor', ones(1,3)-VisRates(i,timestep), 'edgecolor', 'k');
                
            end
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            ylim([0-0.5 1.2+0.5]);
            xlim([0-0.5 1.2+0.5]);
            title('LANDMARKS');
        end
        
    end
    %plot a quiver
    hq = quiver(pos_x(timestep),pos_y(timestep),head_x(timestep)*0.05,head_y(timestep)*0.05,'Linestyle','none');
    
    %now plot annotation to be rat's body
    X = hq.XData;
    Y = hq.YData;
    U = hq.UData;
    V = hq.VData;
    
    ah = annotation('arrow',...
        'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth, 'Color',[0 0 1 0.8]);
    set(ah,'parent',gca);
    set(ah,'position',[X Y lineLength*U lineLength*V]);
    
    hold on
    if(timestep>5)
        line(pos_x(timestep-5:timestep),pos_y(timestep-5:timestep),'Color',[0 0 1 0.3])
    end
    
    if ~visual_scene
        if compartments > 2.5 %tripolar movie
            x = [0,1.2, 0.6, 0];
            y = [0, 0 , height, 0];
            plot(x, y, 'k-', 'LineWidth', 1.0);
            plot(TopBarrier(1,:), TopBarrier(2,:), 'k-', 'LineWidth', 1.0);
            plot(LeftBarrier(1,:), LeftBarrier(2,:), 'k-', 'LineWidth', 1.0);
            plot(RightBarrier(1,:), RightBarrier(2,:), 'k-', 'LineWidth', 1.0);
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            ylim([0-0.5 height+0.5]);
            xlim([0-0.5 1.2+0.5]);
        else  %box-based
            x = [0, 1.2, 1.2, 0];
            y = [0, 0 , 1.2, 1.2];
            plot(x, y, 'k-', 'LineWidth', 1.0);
            plot([0 0],[0 1.2], 'k-','LineWidth',1.0);
            if compartments>1.5
                plot([0.6 0.6],[0 0.55],'k-','LineWidth',1.0);
                plot([0.6 0.6],[0.65 1.2],'k-','LineWidth',1.0);
            end
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            ylim([0-0.5 1.2+0.5]);
            xlim([0-0.5 1.2+0.5]);
        end
    end
    
    % Saving Frame
    suptitle(['Timestep: ',num2str(timestep)]);
    mov(:,timestep)=getframe(1);
    %pause
    pause(0.02);
    
end

movie2avi(mov,'movie.avi','compression','None','quality',100);
end

function rates = file_load(cells, steps, fname)

rates = zeros(cells, steps);

fid = fopen(fname, 'rb');

rates = fread(fid, [steps, cells], 'float32')';

fclose(fid);


end
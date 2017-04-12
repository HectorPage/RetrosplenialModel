function interrogation_movie(seconds,tsize,conjunctive_cells,hd_cells,adn_cells,num_landmarks,start_time, end_time)
%% MOVIE OF NETWORK BEHAVIOUR

disp('Loading Data...')
%% SETTING UP TIMING INFO
steps = uint32(((seconds/tsize)/100)); %Dividing by 100, as new program saves every 100th timestep, to reduce filesize.

if ~start_time %can just set at 0 to go from start
    start_steps = 1;
else
    start_steps = uint32(((start_time/tsize)/100));
end
end_steps = uint32(((end_time/tsize)/100));


vis_cells = num_landmarks; %originally 12 of them


%% Getting rates of the network elements

    ConjunctiveRates = file_load(steps, conjunctive_cells, 'conjunctiveRates_interrogation.bdat')'; %these are transposed, as data output different dimensions to C codes
    HDRates = file_load(steps,hd_cells,'HDRates_interrogation.bdat')';
    VisRates = file_load( steps,vis_cells, 'visRates_interrogation.bdat')';
    ADNRates = file_load(steps,adn_cells,'ADNRates_interrogation.bdat')';

%% Getting directions (PVectors) of each network region

disp('Generating LSE and RSC Population Vectors....')
%Caluclating PVs

%favoured view needed
increment = 360/hd_cells;
favoured_view = (0:hd_cells-1)*increment;
favoured_view = favoured_view';

fav_view_sin = sind(favoured_view); %sine
fav_view_cos = cosd(favoured_view); %cosine

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

vis_panels = [1 2 3 4];
all_conjunctive = [5 6 7 8];
rsc_panels = [9 10 11 12];
ADN_panels = [13 14 15 16];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Visual Cells
figure('doublebuffer','on','units','normalized','outerposition',[0 0 0.5 1]);  %half-screen buffered figure
hold on;
subplot(4,4,vis_panels);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RSCHD Cells
subplot(4,4,rsc_panels);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HD Plot
subplot(4,4,ADN_panels);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Initial Frame
mov=getframe(1);
%% Central Plotting Script

for timestep = start_steps:end_steps
    
    % VIS
    subplot(4,4,vis_panels);
    plot(VisRates(:,timestep),'b','Linewidth',2.0);
    set(gca,'XTick',0:vis_cells);
    hold off
    title('VIS');

    % RSC HD
    subplot(4,4,rsc_panels);
    plot(HDRates(:,timestep),'r','Linewidth',2.0);
    hold on
    plot([RSCPV(timestep)/increment RSCPV(timestep)/increment],[0 1],'Color',[1 0 0],'Linewidth',2.0, 'Linestyle',':');
    plot([ADNPV(timestep)/increment ADNPV(timestep)/increment],[0 1],'Color','m','Linewidth',2.0, 'Linestyle','-');
    set(gca,'XTick',0:hd_cells/4:hd_cells);
    set(gca,'XTickLabel',0:(hd_cells/4)*increment:360);
    hold off
    title('RSC HD');
    
  
        % WC and BC distinction meaningless - plot whole conjunctive layer
        conjunctive_increment = 360/conjunctive_cells;
        subplot(4,4,all_conjunctive);
        plot(ConjunctiveRates(:,timestep),'Color','k','Linewidth',2.0);
        set(gca,'XTick',0:(conjunctive_cells/4):conjunctive_cells);
        set(gca,'XTickLabel',0:(conjunctive_cells/4)*conjunctive_increment:360);
        xlim([0 conjunctive_cells]);
        ylim([0 1.1]);
        hold off
        title('Conjunctive Cells');

    
    % True HD
    subplot(4,4,ADN_panels);
    plot(ADNRates(:,timestep),'m','Linewidth',2.0);
    hold on
   
    plot([RSCPV(timestep)/increment RSCPV(timestep)/increment],[0 1],'Color',[1 0 0],'Linewidth',2.0, 'Linestyle',':');
    plot([ADNPV(timestep)/increment ADNPV(timestep)/increment],[0 1],'Color','m','Linewidth',2.0, 'Linestyle','-');
    set(gca,'XTick',0:hd_cells/4:hd_cells);
    set(gca,'XTickLabel',0:(hd_cells/4)*increment:360);
    hold off
    title(['ADN PV: ', num2str(ADNPV(timestep))]);
    
  
    
    % Saving Frame
    suptitle(['Timestep: ',num2str(timestep)]);
    mov(:,timestep)=getframe(1);
    %pause
    pause(0.02);
    
end

end

function rates = file_load(cells, steps, fname)

rates = zeros(cells, steps);

fid = fopen(fname, 'rb');

rates = fread(fid, [steps, cells], 'float32')';

fclose(fid);


end
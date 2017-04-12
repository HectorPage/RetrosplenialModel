function RSC_model_direction_plot_v3(test_type,full_SO, mode,compartments,seconds, tsize, conjunctive_cells,vis_cells, hd_cells, save_mode,varargin)
%This function shows true Head Direction vs PVector direction of various
%RSC model regions.

%Created by and copyright held by Dr. Hector JI Page 07/11/16
%modified by Dr. Hector JI Page 03/02/2017



%USAGE:
% Make the first argument either 'file' to load in rate data from
%file, or 'simulation' to pass data directly from simulation

%Compartments = 2 for original context box, =3 for triangular
%apparatus, and = 1 for other simulations in one box (vis organisation, and cue
%conflicts)

%If mode = 'file', make useage RSC_model_plot('file', compartments,seconds, tsize, conjunctive_cells,vis_cells,...
%hd_cells, save_mode),(i.e. with varargins being nothing)

%Test type tells the model which sort of test has been done on the data to
%plot:
%1. SO_test = no ADN PI input
%2. disorientation = initial ADN packet in random position
%3. none = no test done

%If mode = 'simluation', make usage RSC_model_plot('file', compartments,seconds, tsize, conjunctive_cells,vis_cells,...
%hd_cells, save_mode, varargins) with varargins as follows:
%   VisRates         = varargin{1};
%   HDRates          = varargin{2};
%   ConjunctiveRates = varargin{3};
%   ADNRates         = varargin{4};

test_type = test_type(1); %only need first letter to discriminate

%% Setting up timesteps
steps = uint32(((seconds/tsize)/100)); %Dividing by 100, as program saves every 100th timestep, to reduce filesize.
xtick = [0,steps/4, steps/2, (steps/4)*3, steps];
xticklabel = [0,seconds/4, seconds/2, (seconds/4)*3, seconds];

%% Reading in rates data
if strcmpi(test_type,'s')
    VisRates = file_load(steps, vis_cells, 'visRates_test.bdat')';
    HDRates = file_load(steps,hd_cells,'HDRates_test.bdat')';
    ConjunctiveRates = file_load(steps,conjunctive_cells, 'conjunctiveRates_test.bdat')';
    ADNRates = file_load(steps,hd_cells,'ADNRates_test.bdat')';
elseif strcmpi(test_type,'d')
    VisRates = file_load(steps, vis_cells, 'visRates_disor.bdat')';
    HDRates = file_load(steps,hd_cells,'HDRates_disor.bdat')';
    ConjunctiveRates = file_load(steps,conjunctive_cells, 'conjunctiveRates_disor.bdat')';
    ADNRates = file_load(steps,hd_cells,'ADNRates_disor.bdat')';
else
    if strcmpi(mode,'file')
        VisRates = file_load(steps, vis_cells, 'visRates.bdat')';
        HDRates = file_load(steps,hd_cells,'HDRates.bdat')';
        ConjunctiveRates = file_load(steps,conjunctive_cells, 'conjunctiveRates.bdat')';
        ADNRates = file_load(steps,hd_cells,'ADNRates.bdat')';
        
    elseif strcmpi(mode, 'simulation')
        VisRates         = varargin{1};
        HDRates          = varargin{2};
        ConjunctiveRates = varargin{3};
        ADNRates         = varargin{4};
        
    else
        error('Mode argument must be ''file'' or ''simulation''!');
    end
end


%% Getting True HD over course of simulation
if strcmpi(test_type,'s')
    fid = fopen('headings_x_test.txt', 'r');
    head_x = fscanf(fid ,'%f\n');
    fclose(fid);
    
    fid = fopen('headings_y_test.txt', 'r');
    head_y = fscanf(fid,'%f\n');
    fclose(fid);
elseif strcmpi(test_type,'d')
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

%normalise all headings
head_lengths = sqrt((head_x.^2)+(head_y.^2));

head_x = head_x./head_lengths;
head_y = head_y./head_lengths;

%Work out true HD over course of simulation
trueHD = atan2d(head_y,head_x);
trueHD(trueHD<0) = trueHD(trueHD<0) + (360);

clear head_x head_y; %no longer need them


%% Caluclating VIS and HD PVectors

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


if ~full_SO %No such thing as vis PV otherwise
    vector_1_vis = sum(bsxfun(@times,VisRates,fav_view_sin),1);
    vector_2_vis = sum(bsxfun(@times,VisRates,fav_view_cos),1);
    
    %Get atan2 (four-quadrant corrected) angles
    VisPV = atan2d(vector_1_vis,vector_2_vis);
    
    %Get in the right range (0-360)
    VisPV(VisPV<0) = VisPV(VisPV<0) + 360;
end

%% Calculating ADN PVector
vector_1_adn = sum(bsxfun(@times,ADNRates,fav_view_sin),1);
vector_2_adn = sum(bsxfun(@times,ADNRates,fav_view_cos),1);
%Get atan2 (four-quadrant corrected) angles
ADNPV = atan2d(vector_1_adn,vector_2_adn);
%Get in the right range (0-360)
ADNPV(ADNPV<0) = ADNPV(ADNPV<0)  + 360;

%% Calculating Conjunctive Cell PVs if appropriate

if compartments>1.5 %only for bipolar and tripolar simulations
    %% Working out PVector of WC and BC flip cells
    fav_view_sin = sind(favoured_view); %sine
    fav_view_cos = cosd(favoured_view); %cosine
    
    WCRates = ConjunctiveRates(1:(conjunctive_cells/2),:);
    BCRates  =ConjunctiveRates((conjunctive_cells/2)+1:end,:);
    
    vector_1_WC = sum(bsxfun(@times,WCRates,fav_view_sin),1);
    vector_2_WC = sum(bsxfun(@times,WCRates,fav_view_cos),1);
    
    vector_1_BC = sum(bsxfun(@times,BCRates,fav_view_sin),1);
    vector_2_BC = sum(bsxfun(@times,BCRates,fav_view_cos),1);
    
    %Get atan2 (four-quadrant corrected) angles
    WCPV = atan2d(vector_1_WC,vector_2_WC);
    BCPV = atan2d(vector_1_BC,vector_2_BC);
    
    %Get in the right range (0-360)
    WCPV(WCPV<0) = WCPV(WCPV<0) + 360;
    BCPV(BCPV<0) = BCPV(BCPV<0)  + 360;
end




%% Working out when the VIS is not firing, and not displaying it
% This is done because looking at what direction RSC indicates when
% landmarks *aren't* visible isn't a fair measure
no_landmark_filter = logical(~sum(VisRates,1));
RSCPV(no_landmark_filter) = NaN;
ADNPV(no_landmark_filter) = NaN;

%% Finding out what compartment the model thinks Rat is in

if(compartments>1.5) %If NOT in a one-compartment simulation
    fid = fopen('odour_context.txt', 'r');
    odour_context = fscanf(fid ,'%f\n');
    fclose(fid);
    
    transition = zeros(numel(odour_context),1);
    for idx = 2:numel(transition)
        if(odour_context(idx)~=odour_context(idx-1))
            transition(idx) = 1;
        end
    end
    
    
    %now find times of transitions
    transition_times = find(transition);
end


%% Plotting the directions indicated by various network areas now
figure()
%1. Actual HD of the animal
difference = abs(diff(trueHD));
jdx = 1;
for idx = 1:numel(trueHD)-1
    if(difference(idx)>180) %detect crossings
        %crossings(jdx) = idx-1; %original version, a typo
        crossings(jdx) = idx;
        jdx = jdx+1;
    end
end

hold on
if exist('crossings','var')
    plot(1:crossings(1),trueHD(1:crossings(1)),'k', 'LineWidth',3.0)
    plot(crossings(end):numel(trueHD),trueHD(crossings(end):numel(trueHD)),'k', 'LineWidth',3.0)
    
    for idx = 1:(numel(crossings)-1)
        plot(crossings(idx)+1:crossings(idx+1),trueHD(crossings(idx)+1:crossings(idx+1)),'k', 'LineWidth',3.0)
    end
    clear crossings;
else
    plot(trueHD,'k','LineWidth',3.0);  %un-wraparounded command
end

if ~full_SO
    % 2. HD as indicated by vision (i.e. flipping between compartments)
    difference = abs(diff(VisPV));
    jdx = 1;
    for idx = 1:numel(VisPV)-1
        if(difference(idx)>180) %detect crossings
            %crossings(jdx) = idx-1; %original version, a typo
            crossings(jdx) = idx;
            jdx = jdx+1;
        end
    end
    
    hold on
    if exist('crossings','var')
        plot(1:crossings(1),VisPV(1:crossings(1)),'r--', 'LineWidth',2.0)
        plot(crossings(end):numel(VisPV),VisPV(crossings(end):numel(VisPV)),'r--', 'LineWidth',2.0)
        
        for idx = 1:(numel(crossings)-1)
            plot(crossings(idx)+1:crossings(idx+1),VisPV(crossings(idx)+1:crossings(idx+1)),'r--', 'LineWidth',2.0)
        end
        clear crossings;
    else
        plot(VisPV,'r--','LineWidth',3.0);  %un-wraparounded command
    end
end
%3. RSC HD direction
difference = abs(diff(RSCPV));
jdx = 1;
for idx = 1:numel(RSCPV)-1
    if(difference(idx)>180) %detect crossings
        %crossings(jdx) = idx-1; %original version, a typo
        crossings(jdx) = idx;
        jdx = jdx+1;
    end
end

hold on
if exist('crossings','var')
    plot(1:crossings(1),RSCPV(1:crossings(1)),'b--', 'LineWidth',2.0)
    plot(crossings(end):numel(RSCPV),RSCPV(crossings(end):numel(RSCPV)),'b--', 'LineWidth',2.0)
    
    for idx = 1:(numel(crossings)-1)
        plot(crossings(idx)+1:crossings(idx+1),RSCPV(crossings(idx)+1:crossings(idx+1)),'b--', 'LineWidth',2.0)
    end
    clear crossings;
else
    plot(RSCPV,'b--','Linewidth',3.0);  %un-wraparounded command
end

if compartments>1.5
    %4. Within compartment directions
    difference = abs(diff(WCPV));
    jdx = 1;
    for idx = 1:numel(WCPV)-1
        if(difference(idx)>180) %detect crossings
            %crossings(jdx) = idx-1; %original version, a typo
            crossings(jdx) = idx;
            jdx = jdx+1;
        end
    end
    
    hold on
    if exist('crossings','var')
        plot(1:crossings(1),WCPV(1:crossings(1)),'g--', 'LineWidth',2.0)
        plot(crossings(end):numel(WCPV),WCPV(crossings(end):numel(WCPV)),'g--', 'LineWidth',2.0)
        
        for idx = 1:(numel(crossings)-1)
            plot(crossings(idx)+1:crossings(idx+1),WCPV(crossings(idx)+1:crossings(idx+1)),'g--', 'LineWidth',2.0)
        end
        clear crossings;
    else
        plot(WCPV,'g--','Linewidth',3.0);  %un-wraparounded command
    end
    
    
    %5. Between compartment
    difference = abs(diff(BCPV));
    jdx = 1;
    for idx = 1:numel(BCPV)-1
        if(difference(idx)>180) %detect crossings
            %crossings(jdx) = idx-1; %original version, a typo
            crossings(jdx) = idx;
            jdx = jdx+1;
        end
    end
    
    hold on
    if exist('crossings','var')
        plot(1:crossings(1),BCPV(1:crossings(1)),'c--', 'LineWidth',2.0,'color','c')
        plot(crossings(end):numel(BCPV),BCPV(crossings(end):numel(BCPV)),'c--', 'LineWidth',2.0)
        
        for idx = 1:(numel(crossings)-1)
            plot(crossings(idx)+1:crossings(idx+1),BCPV(crossings(idx)+1:crossings(idx+1)),'c--', 'LineWidth',2.0)
        end
        clear crossings;
    else
        plot(BCPV,'c--','Linewidth',3.0);  %un-wraparounded command
    end
end

ylim([0 360]);
ylabel('Direction');
set(gca,'YTick',0:60:360);
xlabel('Time (s)');
xlim([0 steps]);
set(gca,'XTick',xtick);
set(gca,'XTickLabel',xticklabel);

%6. Direction as indicated by the ADN layer

difference = abs(diff(ADNPV));
jdx = 1;
for idx = 1:numel(ADNPV)-1
    if(difference(idx)>180) %detect crossings
        %crossings(jdx) = idx-1; %original version, a typo
        crossings(jdx) = idx;
        jdx = jdx+1;
    end
end

hold on
if exist('crossings','var')
    plot(1:crossings(1),ADNPV(1:crossings(1)),'m', 'LineWidth',3.0)
    plot(crossings(end):numel(ADNPV),ADNPV(crossings(end):numel(ADNPV)),'m', 'LineWidth',3.0)
    
    for idx = 1:(numel(crossings)-1)
        plot(crossings(idx)+1:crossings(idx+1),ADNPV(crossings(idx)+1:crossings(idx+1)),'m', 'LineWidth',3.0)
    end
    clear crossings;
else
    plot(ADNPV,'m','LineWidth',3.0);  %un-wraparounded command
end

if(compartments>1.5) %If NOT in a one-compartment simulation
    % 7. Now drawing odour boxes for compartments
    %some simulations, odour doesn't change as a test condition
    if(sum(odour_context)==numel(odour_context)) %if odour is just vanilla
        x = [0 numel(odour_context) numel(odour_context) 0];
        y = [0 0 360 360];
        all_patch = patch(x,y,'w','FaceColor','none');
        set(all_patch,'FaceColor','y','FaceAlpha',0.3);
        
    elseif(sum(odour_context)>1) %if odour is just lemon, no need to do anything
        %only do transitions if there is an actual odour change
        for idx = 1:numel(transition_times)-1
            y = [0 0 360 360];
            x = [transition_times(idx) transition_times(idx+1)-1 transition_times(idx+1)-1 transition_times(idx)];
            p = patch(x,y,'w','FaceColor','none');
            set(p,'EdgeColor','none');
            %set colour
            if(odour_context(transition_times(idx))>1.5) %if apple then set to green
                set(p,'FaceColor','g','FaceAlpha',0.3);
            elseif(odour_context(transition_times(idx)))%if vanilla then set to yellow
                set(p,'FaceColor','y','FaceAlpha',0.3);
            end
        end
        %draw the initial odour
        x_start = [0 transition_times(1)-1 transition_times(1)-1 0];
        y_start = [0 0 360 360];
        p_start = patch(x_start,y_start,'w','FaceColor','none');
        set(p_start,'EdgeColor','none');
        %set colour
        if(odour_context(1)>1.5) %if apple then set to green
            set(p_start,'FaceColor','g','FaceAlpha',0.3);
        elseif(odour_context(1))%if vanilla then set to yellow
            set(p_start,'FaceColor','y','FaceAlpha',0.3);
        end
        
        %draw the final odour
        x_end = [transition_times(end) max(xlim) max(xlim) transition_times(end)];
        y_end = [0 0 360 360];
        p_end = patch(x_end,y_end,'w','FaceColor','none');
        set(p_end,'EdgeColor','none');
        %set colour
        if(odour_context(end)>1.5) %if apple then set to green
            set(p_end,'FaceColor','g','FaceAlpha',0.3);
        elseif(odour_context(end))%if vanilla then set to yellow
            set(p_end,'FaceColor','y','FaceAlpha',0.3);
        end
    end
end

%% Create custom legend (avoids issue of wraparound using multiple plot commands)

if compartments>1.5
    h = zeros(5, 1);
    h(1) = plot(0,0,'k-', 'visible', 'off');
    h(2) = plot(0,0,'r--', 'visible', 'off');
    h(3) = plot(0,0,'b--', 'visible', 'off');
    h(4) = plot(0,0,'g--', 'visible', 'off');
    h(5) = plot(0,0,'c--', 'visible', 'off');
    h(6) = plot(0,0,'m', 'visible', 'off');
    legend(h,{'True HD','VIS Direction','RSC HD', 'WC', 'BC','ADN HD'});
elseif ~full_SO
    h = zeros(3, 1);
    h(1) = plot(0,0,'k-', 'visible', 'off');
    h(2) = plot(0,0,'r--', 'visible', 'off');
    h(3) = plot(0,0,'b--', 'visible', 'off');
    h(4) = plot(0,0,'m', 'visible', 'off');
    legend(h,{'True HD','VIS Direction','RSC HD','ADN HD'});
else
    h = zeros(2, 1);
    h(1) = plot(0,0,'k-', 'visible', 'off');
    h(2) = plot(0,0,'b--', 'visible', 'off');
    h(3) = plot(0,0,'m', 'visible', 'off');
    legend(h,{'HD','RSC HD','ADN HD'});
end


    
if (strcmpi(save_mode, 'save'))
    if strcmpi(test_type,'s')
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        saveas(gcf,'_PVECTORS_test', 'fig');
        close(gcf);
    elseif strcmpi(test_type,'d')
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        saveas(gcf,'_PVECTORS_disor', 'fig');
        close(gcf);
    else
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        saveas(gcf,'_PVECTORS', 'fig');
        close(gcf);
    end
end


    tensecs = uint32(((10.0/tsize)/100)); %Dividing by 100, as program saves every 100th timestep, to reduce filesize.

    RSC_inaccuracy = atan2d(sind(trueHD(1:end-1)-RSCPV'),cosd(trueHD(1:end-1)-RSCPV'));
    %RSC_inaccuracy = abs(RSC_inaccuracy); %want to know about dir
    
     if strcmpi(test_type,'s')
        fid = fopen('_RSC_tracking_stats_test.txt','w');
     elseif strcmpi(test_type,'d')
        fid = fopen('_RSC_tracking_stats_disor.txt','w');
     else
         fid = fopen('_RSC_tracking_stats.txt','w');
     end
    fprintf(fid,'Mean inaccuracy: %f deg\n  SD inaccuracy: %f deg\n In range %f deg to %f deg',nanmean(RSC_inaccuracy(end-tensecs+1:end)),...
        nanstd(RSC_inaccuracy(end-tensecs+1:end)), min(RSC_inaccuracy(end-tensecs+1:end)),max(RSC_inaccuracy(end-tensecs+1:end)));
    fclose(fid);
    
    

    ADN_inaccuracy = atan2d(sind(trueHD(1:end-1)-ADNPV'),cosd(trueHD(1:end-1)-ADNPV'));
    
    
    if strcmpi(test_type,'s')    
           fid = fopen('_ADN_tracking_stats_test.txt','w');
     elseif strcmpi(test_type,'d')
        fid = fopen('_ADN_tracking_stats_disor.txt','w');
     else
         fid = fopen('_ADN_tracking_stats.txt','w');
     end
     fprintf(fid,'Mean inaccuracy: %f deg\n  SD inaccuracy: %f deg\n In range %f deg to %f deg',nanmean(ADN_inaccuracy(end-tensecs+1:end)),...
         nanstd(ADN_inaccuracy(end-tensecs+1:end)),min(ADN_inaccuracy(end-tensecs+1:end)),max(ADN_inaccuracy(end-tensecs+1:end)));
     
     if(sum(ADNRates(:,1)<0.01))
         fprintf(fid,'\n ADN not active');
     end
     
     fclose(fid);
     
    figure()
    plot(RSC_inaccuracy);
    
    if(sum(ADNRates(:,1)>0.01))
        hold on
        plot(ADN_inaccuracy);
    end
    title('RSC/ADN Inaccuracy');
    ylim([-180 180]);
    ylabel('Inaccuracy (deg)');
    xlabel('Timestep');
    xlim([0 steps]);
    
    if (strcmpi(save_mode, 'save'))
        if strcmpi(test_type,'s')  
            set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
            saveas(gcf,'_RSCADN_inaccuracy_test', 'fig');
            close(gcf);
        elseif strcmpi(test_type,'d')  
            set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
            saveas(gcf,'_RSCADN_inaccuracy_disor', 'fig');
            close(gcf);
        else
            set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
            saveas(gcf,'_RSCADN_inaccuracy', 'fig');
            close(gcf);
        end
    end
    

end
function rates = file_load(cells, steps, fname)

rates = zeros(cells, steps);

fid = fopen(fname, 'rb');

rates = fread(fid, [steps, cells], 'float32')';

fclose(fid);


end
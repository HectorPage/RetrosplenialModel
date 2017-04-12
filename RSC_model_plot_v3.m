function RSC_model_plot_v3(test_type,mode,compartments,seconds, tsize, conjunctive_cells,vis_cells, hd_cells, save_mode,varargin)
%Plots lots of simulation info for RSC/conjunctive cell model in various
%environments

%Created by and copyright held by Dr. Hector Page 07/10/16
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

%If mode = 'simulation', make usage RSC_model_plot('file', compartments,seconds, tsize, conjunctive_cells,vis_cells,...
%hd_cells, save_mode, varargins) with varargins as follows:
%   VisRates         = varargin{1};
%   HDRates          = varargin{2};
%   ConjunctiveRates = varargin{3};
%   ADNRates         = varargin{4};

test_type = test_type(1);

%% Determining whether to read in data or have it passed directly from the simulation
steps = uint32(((seconds/tsize)/100)); %Dividing by 100, as program saves every 100th timestep, to reduce filesize.

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
        error('First argument must be ''file'' or ''simulation''!');
    end
end

%% Plot Rates
    %Vis Rates
    plotRates(seconds,tsize,VisRates,'VIS',vis_cells, 'save',test_type);
    %HD Rates
    plotRates(seconds,tsize,HDRates,'HD',hd_cells, 'save',test_type);
    %ADN Rates
    plotRates(seconds,tsize,ADNRates,'ADN',hd_cells, 'save',test_type);
    %Conjunctive Rates
    if(compartments > 1)
        plotRates(seconds,tsize,ConjunctiveRates(1:conjunctive_cells/2,:),'WC',conjunctive_cells/2,'save',test_type);
        plotRates(seconds,tsize,ConjunctiveRates(conjunctive_cells/2+1:end,:),'BC',conjunctive_cells/2,'save',test_type);
    else
        plotRates(seconds,tsize,ConjunctiveRates,'CONJ',conjunctive_cells,'save',1);
    end

%% Plotting PVectors vs True HD
    RSC_model_direction_plot_v3(test_type,1,'simulation',1,seconds,tsize,conjunctive_cells,vis_cells,hd_cells,'save',...
        VisRates,HDRates,ConjunctiveRates,ADNRates);


if compartments > 1.5
    %% Need odour context
    fid = fopen('odour_context.txt', 'r');
    odour_context = fscanf(fid ,'%f\n');
    fclose(fid);
    
    %N.B. Not set to make a difference for SO_test or not: only one
    %compartment simulations so far
end
%% Need trueHD
%% Getting True HD over course of simulation
if strcmpi(test_type,'s')
    fid = fopen('headings_x_test.txt', 'r');
    head_x = fscanf(fid ,'%f\n');
    fclose(fid);
    
    fid = fopen('headings_y_test.txt', 'r');
    head_y = fscanf(fid,'%f\n');
    fclose(fid);
    
    fid = fopen('positions_x_test.txt', 'r');
    pos_x = fscanf(fid ,'%f\n');
    fclose(fid);
    
    fid = fopen('positions_y_test.txt', 'r');
    pos_y = fscanf(fid,'%f\n');
    fclose(fid);
    
    if vis_cells>12 %these are the quadrant filters for a circular apparatus
        TR_filter = pos_x>=0 & pos_y>=0;
        TR_filter(end) = [];

        BR_filter = pos_x>=0 & pos_y<0;
        BR_filter(end) = [];

        TL_filter = pos_x<0 & pos_y>=0;
        TL_filter(end) = [];

        BL_filter = pos_x<0 & pos_y<0;
        BL_filter(end) = [];
    else
        TR_filter = pos_x>=0.6 & pos_y>=0.6;
        TR_filter(end) = [];

        BR_filter = pos_x>=0.6 & pos_y<0.6;
        BR_filter(end) = [];

        TL_filter = pos_x<0.6 & pos_y>=0.6;
        TL_filter(end) = [];

        BL_filter = pos_x<0.6 & pos_y<0.6;
        BL_filter(end) = [];
    end
    
    nBins_shift = 61;
    topEdge_shift = 183;
    bottomEdge_shift = -183;
    binEdges_shift = linspace(bottomEdge_shift, topEdge_shift, nBins_shift+1);
    
elseif  strcmpi(test_type,'d')
    fid = fopen('headings_x_disor.txt', 'r');
    head_x = fscanf(fid ,'%f\n');
    fclose(fid);
    
    fid = fopen('headings_y_disor.txt', 'r');
    head_y = fscanf(fid,'%f\n');
    fclose(fid);
    
    fid = fopen('positions_x_disor.txt', 'r');
    pos_x = fscanf(fid ,'%f\n');
    fclose(fid);
    
    fid = fopen('positions_y_disor.txt', 'r');
    pos_y = fscanf(fid,'%f\n');
    fclose(fid);
    
    if vis_cells>12 %these are the quadrant filters for a circular apparatus
        TR_filter = pos_x>=0 & pos_y>=0;
        TR_filter(end) = [];

        BR_filter = pos_x>=0 & pos_y<0;
        BR_filter(end) = [];

        TL_filter = pos_x<0 & pos_y>=0;
        TL_filter(end) = [];

        BL_filter = pos_x<0 & pos_y<0;
        BL_filter(end) = [];
    else
        TR_filter = pos_x>=0.6 & pos_y>=0.6;
        TR_filter(end) = [];

        BR_filter = pos_x>=0.6 & pos_y<0.6;
        BR_filter(end) = [];

        TL_filter = pos_x<0.6 & pos_y>=0.6;
        TL_filter(end) = [];

        BL_filter = pos_x<0.6 & pos_y<0.6;
        BL_filter(end) = [];
    end
    
    nBins_shift = 61;
    topEdge_shift = 183;
    bottomEdge_shift = -183;
    binEdges_shift = linspace(bottomEdge_shift, topEdge_shift, nBins_shift+1);
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
trueHD = atan2(head_y,head_x);
trueHD(trueHD<0) = trueHD(trueHD<0) + (2*pi);
trueHD = trueHD * (180/pi);

trueHD(end) = []; %it's one element too long

clear head_x head_y; %no longer need them

%% Working out conjunctive cell tuning curves
nBins = 60;
topEdge = 360;
bottomEdge = 0;

binEdges = linspace(bottomEdge, topEdge, nBins+1);

the_angles = 3:6:357;
angles = zeros(nBins+1,1);
angles(1:end-1) = the_angles;
angles(end) = the_angles(1);
angles = angles .* (pi/180);

five_secs = uint32(((5.0/tsize)/100)); %five seconds worth of timesteps

%%N.B.CHECK SO TEST CAPABILITY FOR MULTI COMPARTMENT


if(compartments>1.5) %in two and three compartment simulations, get tuning curve in each compartment
    BC_corrpks_lem = zeros(conjunctive_cells/2,3); %first element is num peaks, second is first peak dist, third is second peak distance
    WC_corrpks_lem = zeros(conjunctive_cells/2,3);
    
    BC_corrpks_van = zeros(conjunctive_cells/2,3);
    WC_corrpks_van = zeros(conjunctive_cells/2,3);
    
    BC_corrpks_all = zeros(conjunctive_cells/2,3);
    WC_corrpks_all = zeros(conjunctive_cells/2,3);
    
    WC_lemon_tuning_curve = zeros(conjunctive_cells/2,61);
    WC_vanilla_tuning_curve = zeros(conjunctive_cells/2,61);
    WC_overall_tuning_curve = zeros(conjunctive_cells/2,61);
    
    BC_lemon_tuning_curve = zeros(conjunctive_cells/2,61);
    BC_vanilla_tuning_curve = zeros(conjunctive_cells/2,61);
    BC_overall_tuning_curve = zeros(conjunctive_cells/2,61);
    
    WC_lemon_autocorr = zeros(conjunctive_cells/2,61);
    WC_vanilla_autocorr = zeros(conjunctive_cells/2,61);
    WC_overall_autocorr = zeros(conjunctive_cells/2,61);
    
    BC_lemon_autocorr = zeros(conjunctive_cells/2,61);
    BC_vanilla_autocorr = zeros(conjunctive_cells/2,61);
    BC_overall_autocorr = zeros(conjunctive_cells/2,61);
    
    if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
        BC_corrpks_app = zeros(conjunctive_cells/2,3);
        WC_corrpks_app = zeros(conjunctive_cells/2,3);
        
        WC_apple_tuning_curve = zeros(conjunctive_cells/2,61);
        BC_apple_tuning_curve = zeros(conjunctive_cells/2,61);
        
        WC_apple_autocorr = zeros(conjunctive_cells/2,61);
        BC_apple_autocorr = zeros(conjunctive_cells/2,61);
        
        apple_condition = odour_context>1.5;
        apple_condition = apple_condition(steps-(five_secs*2)+1:steps);
        trueHD_apple = trueHD(steps-(five_secs*2)+1:steps);
        trueHD_apple = trueHD_apple(apple_condition);
    end
    
    
    %JUST THE LAST 10s OF THE TRIAL
    odour_context = odour_context(steps-(five_secs*2)+1:steps);
    vanilla_condition = odour_context>0.5 & odour_context<1.5;
    lemon_condition = odour_context<0.5;
    
    
    trueHD_vanilla = trueHD(steps-(five_secs*2)+1:steps);
    trueHD_vanilla = trueHD_vanilla(vanilla_condition);
    trueHD_lemon = trueHD(steps-(five_secs*2)+1:steps);
    trueHD_lemon = trueHD_lemon(lemon_condition);
    
    trueHD_all = trueHD(steps-(five_secs*2)+1:steps);
    
    for idx = 1:conjunctive_cells/2
        
        WC_cell = idx;
        BC_cell = idx+(conjunctive_cells/2);
        
        %DOING INDIVIDUAL COMPARTMENTS
        WC_rates_vanilla = ConjunctiveRates(WC_cell,steps-(five_secs*2)+1:steps);
        WC_rates_vanilla = WC_rates_vanilla(vanilla_condition);
        BC_rates_vanilla = ConjunctiveRates(BC_cell,steps-(five_secs*2)+1:steps);
        BC_rates_vanilla = BC_rates_vanilla(vanilla_condition);
        
        
        WC_rates_lemon = ConjunctiveRates(WC_cell,steps-(five_secs*2)+1:steps);
        WC_rates_lemon = WC_rates_lemon(lemon_condition);
        BC_rates_lemon = ConjunctiveRates(BC_cell,steps-(five_secs*2)+1:steps);
        BC_rates_lemon = BC_rates_lemon(lemon_condition);
        
        
        WC_van = generate_tuning_curve(trueHD_vanilla,WC_rates_vanilla,nBins,binEdges);
        WC_vanilla_tuning_curve(idx,:) = WC_van;
        WC_lem = generate_tuning_curve(trueHD_lemon,WC_rates_lemon,nBins,binEdges);
        WC_lemon_tuning_curve(idx,:) = WC_lem;
        
        BC_van = generate_tuning_curve(trueHD_vanilla,BC_rates_vanilla,nBins,binEdges);
        BC_vanilla_tuning_curve(idx,:) = BC_van;
        BC_lem = generate_tuning_curve(trueHD_lemon,BC_rates_lemon,nBins,binEdges);
        BC_lemon_tuning_curve(idx,:) = BC_lem;
        
        if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
            
            WC_rates_apple = ConjunctiveRates(WC_cell,steps-(five_secs*2)+1:steps);
            WC_rates_apple = WC_rates_apple(apple_condition);
            BC_rates_apple = ConjunctiveRates(BC_cell,steps-(five_secs*2)+1:steps);
            BC_rates_apple = BC_rates_apple(apple_condition);
            
            WC_app = generate_tuning_curve(trueHD_apple,WC_rates_apple,nBins,binEdges);
            WC_apple_tuning_curve(idx,:) = WC_app;
            BC_app = generate_tuning_curve(trueHD_apple,BC_rates_apple,nBins,binEdges);
            BC_apple_tuning_curve(idx,:) = BC_app;
        end
        
        
        %ALSO DOING ALL COMPARTMENTS
        WC_rates = ConjunctiveRates(WC_cell,steps-(five_secs*2)+1:steps);
        BC_rates = ConjunctiveRates(BC_cell,steps-(five_secs*2)+1:steps);
        
        WC_all = generate_tuning_curve(trueHD_all,WC_rates,nBins,binEdges);
        WC_overall_tuning_curve(idx,:) = WC_all;
        BC_all = generate_tuning_curve(trueHD_all,BC_rates,nBins,binEdges);
        BC_overall_tuning_curve(idx,:) = BC_all;
        
        %Getting peaks and separations from autocorrelation
        autocorr_threshold = 0.5;
        % Lemon
        
        [BC_lemon_autocorr(idx,:), BC_corrpks_lem(idx,1),BC_corrpks_lem(idx,2:end)] = ...
            tuning_curve_autocorrelation(autocorr_threshold,BC_lem);
        
        
        [WC_lemon_autocorr(idx,:), WC_corrpks_lem(idx,1),WC_corrpks_lem(idx,2:end)] = ...
            tuning_curve_autocorrelation(autocorr_threshold,WC_lem);
        
        
        %Vanilla
        [BC_vanilla_autocorr(idx,:), BC_corrpks_van(idx,1),BC_corrpks_van(idx,2:end)] = ...
            tuning_curve_autocorrelation(autocorr_threshold,BC_van);
        
        
        [WC_vanilla_autocorr(idx,:), WC_corrpks_van(idx,1),WC_corrpks_van(idx,2:end)] = ...
            tuning_curve_autocorrelation(autocorr_threshold,WC_van);
        
        
        %All compartments
        
        [BC_overall_autocorr(idx,:), BC_corrpks_all(idx,1),BC_corrpks_all(idx,2:end)] = ...
            tuning_curve_autocorrelation(autocorr_threshold,BC_all);
        
        
        [WC_overall_autocorr(idx,:), WC_corrpks_all(idx,1),WC_corrpks_all(idx,2:end)] = ...
            tuning_curve_autocorrelation(autocorr_threshold,WC_all);
        %Apple if it exists
        if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
            [BC_apple_autocorr(idx,:), BC_corrpks_app(idx,1),BC_corrpks_app(idx,2:end)] = ...
                tuning_curve_autocorrelation(autocorr_threshold,BC_app);
            
            
            [WC_apple_autocorr(idx,:), WC_corrpks_app(idx,1),WC_corrpks_app(idx,2:end)] = ...
                tuning_curve_autocorrelation(autocorr_threshold,WC_app);
            
        end
        
        if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
            if WC_cell == 150;
                plot_tuning_curve_v2(angles,compartments,odour_context,'save',WC_cell,0,'WC',WC_van,WC_lem,WC_app);
            end
            
            if BC_cell == 350;
                plot_tuning_curve_v2(angles,compartments,odour_context,'save',BC_cell,0,'BC',BC_van,BC_lem,BC_app);
            end
        else
            if WC_cell == 150;
                plot_tuning_curve_v2(angles,compartments,odour_context,'save',WC_cell,0,'WC',WC_van,WC_lem);
            end
            
            if BC_cell == 350;
                plot_tuning_curve_v2(angles,compartments,odour_context,'save',BC_cell,0,'BC',BC_van,BC_lem);
            end
            
        end
    end
    
    
    figure()
    
    cmax = max(max(WC_overall_autocorr));
    
    %PLOTTING AUTOCORRELATIONS
    figure()
    subplot(compartments+1,1,1);
    imagesc(WC_lemon_autocorr)
    set(gca,'YDir','normal')
    set(gca,'YTickLabel',50:50:250);
    set(gca,'XTick',1:10:61);
    set(gca,'XTickLabel',0:60:360);
    set(gca,'TickLength',[ 0 0 ]);
    caxis([0 cmax]);
    set(gca,'Fontsize',14);
    title('WC Lemon');
    
    subplot(compartments+1,1,2);
    imagesc(WC_vanilla_autocorr);
    set(gca,'YDir','normal')
    set(gca,'YTickLabel',50:50:250);
    set(gca,'XTick',1:10:61);
    set(gca,'XTickLabel',0:60:360);
    set(gca,'TickLength',[ 0 0 ]);
    caxis([0 cmax]);
    set(gca,'Fontsize',14);
    title('WC Vanilla');
    
    
    if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
        subplot(compartments+1,1,3);
        imagesc(WC_apple_autocorr);
        set(gca,'YDir','normal')
        set(gca,'YTickLabel',50:50:250);
        set(gca,'XTick',1:10:61);
        set(gca,'XTickLabel',0:60:360);
        set(gca,'TickLength',[ 0 0 ]);
        caxis([0 cmax]);
        set(gca,'Fontsize',14);
        title('WC Apple');
        
        
        subplot(compartments+1,1,4);
        imagesc(WC_overall_autocorr);
        set(gca,'YDir','normal')
        set(gca,'YTickLabel',50:50:250);
        set(gca,'XTick',1:10:61);
        set(gca,'XTickLabel',0:60:360);
        set(gca,'TickLength',[ 0 0 ]);
        caxis([0 cmax]);
        set(gca,'Fontsize',14);
        title('WC All');
        
        
    else
        subplot(compartments+1,1,3);
        imagesc(WC_overall_autocorr);
        set(gca,'YDir','normal')
        set(gca,'YTickLabel',50:50:250);
        set(gca,'XTick',1:10:61);
        set(gca,'XTickLabel',0:60:360);
        set(gca,'TickLength',[ 0 0 ]);
        caxis([0 cmax]);
        set(gca,'Fontsize',14);
        title('WC All');
        
        
    end
    hp4 = get(subplot(compartments+1,1,compartments+1),'Position');
    cbar1 = colorbar('Position', [hp4(1)+hp4(3)+0.02 hp4(2)  0.01  hp4(2)+hp4(3)*0.9]);
    set(cbar1,'Fontsize',14);
    if (strcmpi(save_mode, 'save'))
        if strcmpi(test_type,'s')
            saveas(gcf,'_WC_autocorr_test','fig');
        elseif strcmpi(test_type,'d')
            saveas(gcf,'_WC_autocorr_disor','fig');
        else
            saveas(gcf,'_WC_autocorr', 'fig');
            close(gcf);
        end
        
    end
    
    cmax = max(max(BC_overall_autocorr));
    
    figure()
    subplot(compartments+1,1,1);
    imagesc(BC_lemon_autocorr);
    set(gca,'YDir','normal')
    set(gca,'YTickLabel',50:50:250);
    set(gca,'XTick',1:10:61);
    set(gca,'XTickLabel',0:60:360);
    set(gca,'TickLength',[ 0 0 ]);
    caxis([0 cmax]);
    title('BC Lemon');
    set(gca,'Fontsize',14);
    
    
    subplot(compartments+1,1,2);
    imagesc(BC_vanilla_autocorr);
    set(gca,'YDir','normal')
    set(gca,'YTickLabel',50:50:250);
    set(gca,'XTick',1:10:61);
    set(gca,'XTickLabel',0:60:360);
    set(gca,'TickLength',[ 0 0 ]);
    caxis([0 cmax]);
    set(gca,'Fontsize',14);
    title('BC Vanilla');
    
    
    if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
        subplot(compartments+1,1,3);
        imagesc(BC_apple_autocorr);
        set(gca,'YDir','normal')
        set(gca,'YTickLabel',50:50:250);
        set(gca,'XTick',1:10:61);
        set(gca,'XTickLabel',0:60:360);
        set(gca,'TickLength',[ 0 0 ]);
        caxis([0 cmax]);
        set(gca,'Fontsize',14);
        title('BC Apple');
        
        subplot(compartments+1,1,4);
        imagesc(BC_overall_autocorr);
        set(gca,'YDir','normal')
        set(gca,'YTickLabel',50:50:250);
        set(gca,'XTick',1:10:61);
        set(gca,'XTickLabel',0:60:360);
        set(gca,'TickLength',[ 0 0 ]);
        caxis([0 cmax]);
        set(gca,'Fontsize',14);
        title('BC All');
        
    else
        subplot(compartments+1,1,3);
        imagesc(BC_overall_autocorr);
        set(gca,'YDir','normal')
        set(gca,'YTickLabel',50:50:250);
        set(gca,'XTick',1:10:61);
        set(gca,'XTickLabel',0:60:360);
        set(gca,'TickLength',[ 0 0 ]);
        caxis([0 cmax]);
        set(gca,'Fontsize',14);
        title('BC All');
        
        
    end
    hp4 = get(subplot(compartments+1,1,compartments+1),'Position');
    cbar2 = colorbar('Position', [hp4(1)+hp4(3)+0.02 hp4(2)  0.01  hp4(2)+hp4(3)*0.9]);
    set(cbar2,'Fontsize',14);
    
   if (strcmpi(save_mode, 'save'))
        if strcmpi(test_type,'s')
            saveas(gcf,'_BC_autocorr_test','fig');
        elseif strcmpi(test_type,'d')
            saveas(gcf,'_BC_autocorr_disor','fig');
        else
            saveas(gcf,'_BC_autocorr', 'fig');
            close(gcf);
        end
        
    end
    
    
    
    %% Plotting Peak Counts and separations
    %BC
    figure()
    subplot((compartments+1),2,1)
    plot(BC_corrpks_lem(:,1),'Linewidth',2.0);
    xlabel('Cell');
    ylabel('Num. Peaks')
    title('BC Lemon Peaks');
    set(gca,'Fontsize',14);
    ylim([0 3]);
    subplot((compartments+1),2,2)
    plot(BC_corrpks_lem(:,2),'Linewidth',2.0);
    hold on
    plot(BC_corrpks_lem(:,3),'Linewidth',2.0);
    title('Lemon Peak Separation');
    xlabel('Cell');
    ylabel('Peak Separation')
    set(gca,'Fontsize',14);
    ylim([0 200]);
    subplot((compartments+1),2,3)
    plot(BC_corrpks_van(:,1),'Linewidth',2.0);
    title('BC Vanilla Peaks');
    xlabel('Cell');
    ylabel('Num. Peaks')
    set(gca,'Fontsize',14);
    ylim([0 3]);
    subplot((compartments+1),2,4)
    plot(BC_corrpks_van(:,2),'Linewidth',2.0);
    hold on
    plot(BC_corrpks_van(:,3),'Linewidth',2.0);
    title('Vanilla Peak Separation');
    xlabel('Cell');
    ylabel('Peak Separation')
    set(gca,'Fontsize',14);
    ylim([0 200]);
    
    
    if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
        subplot((compartments+1),2,5)
        plot(BC_corrpks_app(:,1),'Linewidth',2.0);
        title('BC Apple Peaks');
        xlabel('Cell');
        ylabel('Num. Peaks')
        set(gca,'Fontsize',14);
        ylim([0 3]);
        subplot((compartments+1),2,6)
        plot(BC_corrpks_app(:,2),'Linewidth',2.0);
        hold on
        plot(BC_corrpks_app(:,3),'Linewidth',2.0);
        title('Apple Peak Separation');
        xlabel('Cell');
        ylabel('Peak Separation')
        set(gca,'Fontsize',14);
        ylim([0 200]);
        subplot((compartments+1),2,7)
        plot(BC_corrpks_all(:,1),'Linewidth',2.0);
        title('BC Overall Peaks');
        xlabel('Cell');
        ylabel('Num. Peaks')
        set(gca,'Fontsize',14);
        ylim([0 3]);
        subplot((compartments+1),2,8)
        plot(BC_corrpks_all(:,2),'Linewidth',2.0);
        hold on
        plot(BC_corrpks_all(:,3),'Linewidth',2.0);
        title('Overall Peak Separation');
        xlabel('Cell');
        ylabel('Peak Separation')
        set(gca,'Fontsize',14);
        ylim([0 200]);
    else
        subplot((compartments+1),2,5)
        plot(BC_corrpks_all(:,1),'Linewidth',2.0);
        title('BC Overall Peaks');
        xlabel('Cell');
        ylabel('Num. Peaks')
        set(gca,'Fontsize',14);
        ylim([0 3]);
        subplot((compartments+1),2,6)
        plot(BC_corrpks_all(:,2),'Linewidth',2.0);
        hold on
        plot(BC_corrpks_all(:,3),'Linewidth',2.0);
        title('Overall Peak Separation');
        xlabel('Cell');
        ylabel('Peak Separation')
        set(gca,'Fontsize',14);
        ylim([0 200]);
    end
    if (strcmpi(save_mode, 'save'))
        if strcmpi(test_type,'s')
            saveas(gcf,'_BC_peaks_test','fig');
        elseif strcmpi(test_type,'d')
            saveas(gcf,'_BC_peaks_disor','fig');
        else
            saveas(gcf,'_BC_peaks', 'fig');
            close(gcf);
        end
    end
    
    %WC
    figure()
    subplot((compartments+1),2,1)
    plot(WC_corrpks_lem(:,1),'Linewidth',2.0);
    xlabel('Cell');
    ylabel('Num. Peaks')
    title('WC Lemon Peaks');
    set(gca,'Fontsize',14);
    ylim([0 3]);
    subplot((compartments+1),2,2)
    plot(WC_corrpks_lem(:,2),'Linewidth',2.0);
    hold on
    plot(WC_corrpks_lem(:,3),'Linewidth',2.0);
    title('Lemon Peak Separation');
    xlabel('Cell');
    ylabel('Peak Separation')
    set(gca,'Fontsize',14);
    ylim([0 200]);
    subplot((compartments+1),2,3)
    plot(WC_corrpks_van(:,1),'Linewidth',2.0);
    title('WC Vanilla Peaks');
    xlabel('Cell');
    ylabel('Num. Peaks')
    set(gca,'Fontsize',14);
    ylim([0 3]);
    subplot((compartments+1),2,4)
    plot(WC_corrpks_van(:,2),'Linewidth',2.0);
    hold on
    plot(WC_corrpks_van(:,3),'Linewidth',2.0);
    title('Vanilla Peak Separation');
    xlabel('Cell');
    ylabel('Peak Separation')
    set(gca,'Fontsize',14);
    ylim([0 200]);
    
    
    if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
        subplot((compartments+1),2,5)
        plot(WC_corrpks_app(:,1),'Linewidth',2.0);
        title('WC Apple Peaks');
        xlabel('Cell');
        ylabel('Num. Peaks')
        ylim([0 3]);
        set(gca,'Fontsize',14);
        subplot((compartments+1),2,6)
        plot(WC_corrpks_app(:,2),'Linewidth',2.0);
        hold on
        plot(WC_corrpks_app(:,3),'Linewidth',2.0);
        title('Apple Peak Separation');
        xlabel('Cell');
        ylabel('Peak Separation')
        set(gca,'Fontsize',14);
        ylim([0 200]);
        subplot((compartments+1),2,7)
        plot(WC_corrpks_all(:,1),'Linewidth',2.0);
        title('WC Overall Peaks');
        xlabel('Cell');
        ylabel('Num. Peaks')
        set(gca,'Fontsize',14);
        ylim([0 3]);
        subplot((compartments+1),2,8)
        plot(WC_corrpks_all(:,2),'Linewidth',2.0);
        hold on
        plot(WC_corrpks_all(:,3),'Linewidth',2.0);
        title('Overall Peak Separation');
        xlabel('Cell');
        ylabel('Peak Separation')
        set(gca,'Fontsize',14);
        ylim([0 200]);
    else
        subplot((compartments+1),2,5)
        plot(WC_corrpks_all(:,1),'Linewidth',2.0);
        title('WC Overall Peaks');
        xlabel('Cell');
        ylabel('Num. Peaks')
        set(gca,'Fontsize',14);
        ylim([0 3]);
        subplot((compartments+1),2,6)
        plot(WC_corrpks_all(:,2),'Linewidth',2.0);
        hold on
        plot(WC_corrpks_all(:,3),'Linewidth',2.0);
        title('Overall Peak Separation');
        xlabel('Cell');
        ylabel('Peak Separation')
        set(gca,'Fontsize',14);
        ylim([0 200]);
    end
    
    if (strcmpi(save_mode, 'save'))
        if strcmpi(test_type,'s')
            saveas(gcf,'_WC_peaks_test','fig');
        elseif  strcmpi(test_type,'d')
            saveas(gcf,'_WC_peaks_disor','fig');
        else
            saveas(gcf,'_WC_peaks', 'fig');
            close(gcf);
        end
    end
    
    %% Plotting the autocorrelation summaries
    figure()
    %WC Lemon
    subplot(compartments+1,1,1);
    [~,~] = boundedline(angles(1:60),mean(WC_lemon_autocorr(:,1:60)),std(WC_lemon_autocorr(:,1:60)));
    xlim([0 2*pi]);
    set(gca,'XTick',0:(2/6)*pi:2*pi);
    set(gca,'XTickLabel',0:60:360);
    ylim([0 12])
    set(gca,'Fontsize',14);
    title('WC Lemon Average');
    xlabel('Rotation');
    ylabel('Correlation');
    
    %WC Vanilla
    subplot(compartments+1,1,2);
    [~,~] = boundedline(angles(1:60),mean(WC_vanilla_autocorr(:,1:60)),std(WC_vanilla_autocorr(:,1:60)));
    xlim([0 2*pi]);
    set(gca,'XTick',0:(2/6)*pi:2*pi);
    set(gca,'XTickLabel',0:60:360);
    ylim([0 12])
    set(gca,'Fontsize',14);
    title('WC Vanilla Average');
    xlabel('Rotation');
    ylabel('Correlation');
    
    if(any(odour_context>1.5))
        %WC Apple
        subplot(compartments+1,1,3);
        [~,~] = boundedline(angles(1:60),mean(WC_apple_autocorr(:,1:60)),std(WC_apple_autocorr(:,1:60)));
        xlim([0 2*pi]);
        set(gca,'XTick',0:(2/6)*pi:2*pi);
        set(gca,'XTickLabel',0:60:360);
        ylim([0 12])
        set(gca,'Fontsize',14);
        title('WC Apple Average');
        xlabel('Rotation');
        ylabel('Correlation');
        
        %WC All
        subplot(compartments+1,1,4);
        [~,~] = boundedline(angles(1:60),mean(WC_overall_autocorr(:,1:60)),std(WC_overall_autocorr(:,1:60)));
        xlim([0 2*pi]);
        set(gca,'XTick',0:(2/6)*pi:2*pi);
        set(gca,'XTickLabel',0:60:360);
        ylim([0 12])
        set(gca,'Fontsize',14);
        title('WC Overall Average');
        xlabel('Rotation');
        ylabel('Correlation');
    else
        %WC All
        subplot(compartments+1,1,3);
        [~,~] = boundedline(angles(1:60),mean(WC_overall_autocorr(:,1:60)),std(WC_overall_autocorr(:,1:60)));
        xlim([0 2*pi]);
        set(gca,'XTick',0:(2/6)*pi:2*pi);
        set(gca,'XTickLabel',0:60:360);
        ylim([0 12])
        set(gca,'Fontsize',14);
        title('WC Overall Average');
        xlabel('Rotation');
        ylabel('Correlation');
    end
    
    
    if (strcmpi(save_mode, 'save'))
        if strcmpi(test_type,'s')
            saveas(gcf,'_WC_autocorr_summary_test', 'fig');
            close(gcf);
        elseif strcmpi(test_type,'d')
            saveas(gcf,'_WC_autocorr_summary_test', 'fig');
            close(gcf);
        else
            saveas(gcf,'_WC_autocorr_summary', 'fig');
            close(gcf);
        end
    end
    
    figure()
    %BC Lemon
    subplot(compartments+1,1,1);
    [~,~] = boundedline(angles(1:60),mean(BC_lemon_autocorr(:,1:60)),std(BC_lemon_autocorr(:,1:60)));
    xlim([0 2*pi]);
    set(gca,'XTick',0:(2/6)*pi:2*pi);
    set(gca,'XTickLabel',0:60:360);
    ylim([0 12])
    set(gca,'Fontsize',14);
    title('BC Lemon Average');
    xlabel('Rotation');
    ylabel('Correlation');
    
    %BC Vanilla
    subplot(compartments+1,1,2);
    [~,~] = boundedline(angles(1:60),mean(BC_vanilla_autocorr(:,1:60)),std(BC_vanilla_autocorr(:,1:60)));
    xlim([0 2*pi]);
    set(gca,'XTick',0:(2/6)*pi:2*pi);
    set(gca,'XTickLabel',0:60:360);
    ylim([0 12])
    set(gca,'Fontsize',14);
    title('BC Vanilla Average');
    xlabel('Rotation');
    ylabel('Correlation');
    
    if(any(odour_context>1.5))
        %BC Apple
        subplot(compartments+1,1,3);
        [~,~] = boundedline(angles(1:60),mean(BC_apple_autocorr(:,1:60)),std(BC_apple_autocorr(:,1:60)));
        xlim([0 2*pi]);
        set(gca,'XTick',0:(2/6)*pi:2*pi);
        set(gca,'XTickLabel',0:60:360);
        ylim([0 12])
        set(gca,'Fontsize',14);
        title('BC Apple Average');
        xlabel('Rotation');
        ylabel('Correlation');
        
        %BC All
        subplot(compartments+1,1,4);
        [~,~] = boundedline(angles(1:60),mean(BC_overall_autocorr(:,1:60)),std(BC_overall_autocorr(:,1:60)));
        xlim([0 2*pi]);
        set(gca,'XTick',0:(2/6)*pi:2*pi);
        set(gca,'XTickLabel',0:60:360);
        ylim([0 12])
        set(gca,'Fontsize',14);
        title('BC Overall Average');
        xlabel('Rotation');
        ylabel('Correlation');
        
    else
        %BC All
        subplot(compartments+1,1,3);
        [~,~] = boundedline(angles(1:60),mean(BC_overall_autocorr(:,1:60)),std(BC_overall_autocorr(:,1:60)));
        xlim([0 2*pi]);
        set(gca,'XTick',0:(2/6)*pi:2*pi);
        set(gca,'XTickLabel',0:60:360);
        ylim([0 12])
        set(gca,'Fontsize',14);
        title('BC Overall Average');
        xlabel('Rotation');
        ylabel('Correlation');
        
    end
    
    if (strcmpi(save_mode, 'save'))
        if strcmpi(test_type,'s')
            saveas(gcf,'_BC_autocorr_summary_test', 'fig');
            close(gcf);
        elseif strcmpi(test_type,'d')
            saveas(gcf,'_BC_autocorr_summary_disor', 'fig');
            close(gcf);
        else
            saveas(gcf,'_BC_autocorr_summary', 'fig');
            close(gcf);
        end
    end
    
    
    if strcmpi(test_type,'s')
        
        %SAVE THE AUTOCORRELATIONS
        save('WC_lemon_autocorr_test.mat','WC_lemon_autocorr');
        save('WC_vanilla_autocorr_test.mat','WC_vanilla_autocorr');
        save('WC_overall_autocorr_test.mat','WC_overall_autocorr');
        
        save('BC_lemon_autocorr_test.mat','BC_lemon_autocorr');
        save('BC_vanilla_autocorr_test.mat','BC_vanilla_autocorr');
        save('BC_overall_autocorr_test.mat','BC_overall_autocorr');
        
        
        if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
            save('WC_apple_autocorr_test.mat','WC_apple_autocorr');
            save('BC_apple_autocorr_test.mat','BC_apple_autocorr');
        end
        
        %SAVE THE CORRELATION STATS
        save('BC_corrpks_lem_test.mat','BC_corrpks_lem');
        save('WC_corrpks_lem_test.mat','WC_corrpks_lem');
        
        save('BC_corrpks_van_test.mat','BC_corrpks_van');
        save('WC_corr_van_test.mat','WC_corrpks_van');
        
        save('BC_corrpks_all_test.mat','BC_corrpks_all');
        save('WC_corrpks_all_test.mat','WC_corrpks_all');
        
        
        if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
            save('BC_corrpks_app_test.mat','BC_corrpks_app');
            save('WC_corrpks_app_test.mat','WC_corrpks_app');
        end
        
        
        %SAVE THE TUNING CURVES
        save('WC_lemon_tuning_curve_test.mat','WC_lemon_tuning_curve');
        save('WC_vanilla_tuning_curve_test.mat','WC_vanilla_tuning_curve');
        save('WC_overall_tuning_curve_test.mat','WC_overall_tuning_curve');
        
        save('BC_lemon_tuning_curve_test.mat','BC_lemon_tuning_curve');
        save('BC_vanilla_tuning_curve_test.mat','BC_vanilla_tuning_curve');
        save('BC_overall_tuning_curve_test.mat','BC_overall_tuning_curve');
        
        if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
            save('WC_apple_tuning_curve_test.mat','WC_apple_tuning_curve');
            save('BC_apple_tuning_curve_test.mat','BC_apple_tuning_curve');
        end
        
    elseif strcmpi(test_type,'d')
        
        %SAVE THE AUTOCORRELATIONS
        save('WC_lemon_autocorr_disor.mat','WC_lemon_autocorr');
        save('WC_vanilla_autocorr_disor.mat','WC_vanilla_autocorr');
        save('WC_overall_autocorr_disor.mat','WC_overall_autocorr');
        
        save('BC_lemon_autocorr_disor.mat','BC_lemon_autocorr');
        save('BC_vanilla_autocorr_disor.mat','BC_vanilla_autocorr');
        save('BC_overall_autocorr_disor.mat','BC_overall_autocorr');
        
        
        if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
            save('WC_apple_autocorr_disor.mat','WC_apple_autocorr');
            save('BC_apple_autocorr_disor.mat','BC_apple_autocorr');
        end
        
        %SAVE THE CORRELATION STATS
        save('BC_corrpks_lem_disor.mat','BC_corrpks_lem');
        save('WC_corrpks_lem_disor.mat','WC_corrpks_lem');
        
        save('BC_corrpks_van_disor.mat','BC_corrpks_van');
        save('WC_corr_van_disor.mat','WC_corrpks_van');
        
        save('BC_corrpks_all_disor.mat','BC_corrpks_all');
        save('WC_corrpks_all_disor.mat','WC_corrpks_all');
        
        
        if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
            save('BC_corrpks_app_disor.mat','BC_corrpks_app');
            save('WC_corrpks_app_disor.mat','WC_corrpks_app');
        end
        
        
        %SAVE THE TUNING CURVES
        save('WC_lemon_tuning_curve_disor.mat','WC_lemon_tuning_curve');
        save('WC_vanilla_tuning_curve_disor.mat','WC_vanilla_tuning_curve');
        save('WC_overall_tuning_curve_disor.mat','WC_overall_tuning_curve');
        
        save('BC_lemon_tuning_curve_disor.mat','BC_lemon_tuning_curve');
        save('BC_vanilla_tuning_curve_disor.mat','BC_vanilla_tuning_curve');
        save('BC_overall_tuning_curve_disor.mat','BC_overall_tuning_curve');
        
        if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
            save('WC_apple_tuning_curve_disor.mat','WC_apple_tuning_curve');
            save('BC_apple_tuning_curve_disor.mat','BC_apple_tuning_curve');
        end
        
    else
        %SAVE THE AUTOCORRELATIONS
        save('WC_lemon_autocorr.mat','WC_lemon_autocorr');
        save('WC_vanilla_autocorr.mat','WC_vanilla_autocorr');
        save('WC_overall_autocorr.mat','WC_overall_autocorr');
        
        save('BC_lemon_autocorr.mat','BC_lemon_autocorr');
        save('BC_vanilla_autocorr.mat','BC_vanilla_autocorr');
        save('BC_overall_autocorr.mat','BC_overall_autocorr');
        
        
        if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
            save('WC_apple_autocorr.mat','WC_apple_autocorr');
            save('BC_apple_autocorr.mat','BC_apple_autocorr');
        end
        
        %SAVE THE CORRELATION STATS
        save('BC_corrpks_lem.mat','BC_corrpks_lem');
        save('WC_corrpks_lem.mat','WC_corrpks_lem');
        
        save('BC_corrpks_van.mat','BC_corrpks_van');
        save('WC_corr_van.mat','WC_corrpks_van');
        
        save('BC_corrpks_all.mat','BC_corrpks_all');
        save('WC_corrpks_all.mat','WC_corrpks_all');
        
        
        if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
            save('BC_corrpks_app.mat','BC_corrpks_app');
            save('WC_corrpks_app.mat','WC_corrpks_app');
        end
        
        
        %SAVE THE TUNING CURVES
        save('WC_lemon_tuning_curve.mat','WC_lemon_tuning_curve');
        save('WC_vanilla_tuning_curve.mat','WC_vanilla_tuning_curve');
        save('WC_overall_tuning_curve.mat','WC_overall_tuning_curve');
        
        save('BC_lemon_tuning_curve.mat','BC_lemon_tuning_curve');
        save('BC_vanilla_tuning_curve.mat','BC_vanilla_tuning_curve');
        save('BC_overall_tuning_curve.mat','BC_overall_tuning_curve');
        
        if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
            save('WC_apple_tuning_curve.mat','WC_apple_tuning_curve');
            save('BC_apple_tuning_curve.mat','BC_apple_tuning_curve');
        end
        
    end
    
else %if only a single compartment simulation
    
    % SET UP THE DATA STORAGE
    
    trueHD_all = trueHD(steps-(five_secs*2)+1:steps);
    corrpks_all = zeros(conjunctive_cells,3);
    overall_autocorr = zeros(conjunctive_cells,61);
    overall_tuning_curve = zeros(conjunctive_cells,61);
    
    
    for conj_cell = 1:conjunctive_cells
        Rates = ConjunctiveRates(conj_cell,steps-(five_secs*2)+1:steps);
        tuning_curve = generate_tuning_curve(trueHD_all,Rates,nBins,binEdges);
        overall_tuning_curve(conj_cell,:) = tuning_curve;
        
        %Getting peaks and separations from autocorrelation
        autocorr_threshold = 0.5;
        
        [overall_autocorr(conj_cell,:), corrpks_all(conj_cell,1),corrpks_all(conj_cell,2:end)] = ...
            tuning_curve_autocorrelation(autocorr_threshold,tuning_curve);
 
            if conj_cell == 150;
                plot_tuning_curve_v2(angles,compartments,0,'save',conj_cell,test_type,'CONJ',tuning_curve);
            end
            
            if conj_cell == 350;
                plot_tuning_curve_v2(angles,compartments,0,'save',conj_cell,test_type,'CONJ',tuning_curve);
            end
        
    end
    
    %Look at average autocorrelations for all CONJ cells
    figure()
    [~,~] = boundedline(angles(1:60),mean(overall_autocorr(:,1:60)),std(overall_autocorr(:,1:60)));
    xlim([0 2*pi]);
    set(gca,'XTick',0:(2/6)*pi:2*pi);
    set(gca,'XTickLabel',0:60:360);
    ylim([0 12])
    set(gca,'Fontsize',14);
    title('All Conj Average');
    xlabel('Rotation');
    ylabel('Correlation');
    if (strcmpi(save_mode, 'save'))
       if strcmpi(test_type,'s')
            saveas(gcf,'_CONJ_autocorr_summary_test', 'fig');
            close(gcf);
        elseif strcmpi(test_type,'d')
            saveas(gcf,'_CONJ_autocorr_summary_disor', 'fig');
            close(gcf);
        else
            saveas(gcf,'_CONJ_autocorr_summary', 'fig');
            close(gcf);
        end
    end
    
    %Look at autocorrelation heat map for all CONJ cells
    figure()
    cmax = max(max(overall_autocorr));
    imagesc(overall_autocorr);
    set(gca,'YDir','normal')
    set(gca,'YTickLabel',conjunctive_cells/10:conjunctive_cells/10:conjunctive_cells);
    set(gca,'XTick',1:10:61);
    set(gca,'XTickLabel',0:60:360);
    set(gca,'TickLength',[ 0 0 ]);
    caxis([0 cmax]);
    set(gca,'Fontsize',14);
    title('All Conj');
    if (strcmpi(save_mode, 'save'))
        if strcmpi(test_type,'s')
            saveas(gcf,'_CONJ_autocorr_test', 'fig');
            close(gcf);
        elseif strcmpi(test_type,'d')
            saveas(gcf,'_CONJ_autocorr_disor', 'fig');
            close(gcf);
        else
            saveas(gcf,'_CONJ_autocorr', 'fig');
            close(gcf);
        end
    end
    
    if  strcmpi(test_type,'s')
        %Saving the data
        save('corrpks_all_test.mat','corrpks_all');
        save('overall_autocorr_test.mat','overall_autocorr');
        save('overall_tuning_curve_test.mat','overall_tuning_curve');
    elseif strcmpi(test_type,'d')
        %Saving the data
        save('corrpks_all_disor.mat','corrpks_all');
        save('overall_autocorr_disor.mat','overall_autocorr');
        save('overall_tuning_curve_disor.mat','overall_tuning_curve');
    else
        %Saving the data
        save('corrpks_all.mat','corrpks_all');
        save('overall_autocorr.mat','overall_autocorr');
        save('overall_tuning_curve.mat','overall_tuning_curve');
    end
    
end

%% NOW DOING OVERALL PLOTS FOR RSC HD AND ADN HD CELLS

if strcmpi(test_type,'s')|| strcmpi(test_type,'d')
    
    %Remember this is only the last 10 seconds of simulation, so edit
    %accordingly
    
    TR_filter = TR_filter(steps-(five_secs*2)+1:steps);
    BR_filter = BR_filter(steps-(five_secs*2)+1:steps);
    TL_filter = TL_filter(steps-(five_secs*2)+1:steps);
    BL_filter = BL_filter(steps-(five_secs*2)+1:steps);
    
    trueHD = trueHD(steps-(five_secs*2)+1:steps);
    
    %Top right quadrant
    trueHD_TR = trueHD(TR_filter);
    shifted_tuning_curves_TR = zeros(hd_cells,61);
    shifted_tuning_curves_TR_ADN = zeros(hd_cells,61);
    
    %Bottom right quadrant
    trueHD_BR = trueHD(BR_filter);
    shifted_tuning_curves_BR = zeros(hd_cells,61);
    shifted_tuning_curves_BR_ADN = zeros(hd_cells,61);
    
    %Top left quadrant
    trueHD_TL = trueHD(TL_filter);
    shifted_tuning_curves_TL = zeros(hd_cells,61);
    shifted_tuning_curves_TL_ADN = zeros(hd_cells,61);
    
    %Bottom left quadrant
    trueHD_BL = trueHD(BL_filter);
    shifted_tuning_curves_BL = zeros(hd_cells,61);
    shifted_tuning_curves_BL_ADN = zeros(hd_cells,61);
    
    [HD_PFDs] = set_PFDs(hd_cells);
    
end

overall_tuning_curve_rschd = zeros(hd_cells,61);
overall_autocorr_rschd = zeros(hd_cells,61);
corrpks_all_rschd = zeros(hd_cells,3);

overall_tuning_curve_adn = zeros(hd_cells,61);
overall_autocorr_adn = zeros(hd_cells,61);
corrpks_all_adn = zeros(hd_cells,3);


for hd_cell = 1:hd_cells
    Rates_HD = HDRates(hd_cell,steps-(five_secs*2)+1:steps);
    Rates_ADN = ADNRates(hd_cell,steps-(five_secs*2)+1:steps);
    
    hd_tuning_curve = generate_tuning_curve(trueHD_all,Rates_HD,nBins,binEdges);
    adn_tuning_curve = generate_tuning_curve(trueHD_all,Rates_ADN,nBins,binEdges);
    
    overall_tuning_curve_rschd(hd_cell,:) = hd_tuning_curve;
    overall_tuning_curve_adn(hd_cell,:) = adn_tuning_curve;
    
    if strcmpi(test_type,'s')|| strcmpi(test_type,'d')
        
        %Top Right Quadrant
            %RSC
        Rates_HD_TR = Rates_HD(TR_filter);
        shifted_HD_TR = trueHD_TR - HD_PFDs(hd_cell);
        shifted_HD_TR(shifted_HD_TR>180) = shifted_HD_TR(shifted_HD_TR>180) - 360;
        hd_tuning_curve_TR = generate_tuning_curve(shifted_HD_TR,Rates_HD_TR,nBins_shift,binEdges_shift);
        hd_tuning_curve_TR(end) = [];
        if any(hd_tuning_curve_TR<0)
            disp('Something odd TR');
        end
        shifted_tuning_curves_TR(hd_cell,:) = hd_tuning_curve_TR;
        
            %ADN
        Rates_ADN_TR = Rates_ADN(TR_filter);
        adn_tuning_curve_TR = generate_tuning_curve(shifted_HD_TR,Rates_ADN_TR,nBins_shift,binEdges_shift);
        adn_tuning_curve_TR(end) = [];
        if any(adn_tuning_curve_TR<0)
            disp('Something odd TR');
        end
        shifted_tuning_curves_TR_ADN(hd_cell,:) = adn_tuning_curve_TR;
        
        %Bottom Right Quadrant
            %RSC
        Rates_HD_BR = Rates_HD(BR_filter);
        shifted_HD_BR = trueHD_BR - HD_PFDs(hd_cell);
        shifted_HD_BR(shifted_HD_BR>180) = shifted_HD_BR(shifted_HD_BR>180) - 360;
        hd_tuning_curve_BR = generate_tuning_curve(shifted_HD_BR,Rates_HD_BR,nBins_shift,binEdges_shift);
        hd_tuning_curve_BR(end) = [];
        if any(hd_tuning_curve_BR<0)
            disp('Something odd BR');
        end
        shifted_tuning_curves_BR(hd_cell,:) = hd_tuning_curve_BR;
        
              %ADN
        Rates_ADN_BR = Rates_ADN(BR_filter);
        adn_tuning_curve_BR = generate_tuning_curve(shifted_HD_BR,Rates_ADN_BR,nBins_shift,binEdges_shift);
        adn_tuning_curve_BR(end) = [];
        if any(adn_tuning_curve_BR<0)
            disp('Something odd BR');
        end
        shifted_tuning_curves_BR_ADN(hd_cell,:) = adn_tuning_curve_BR;
        
        %Top Left Quadrant
            %RSC
        Rates_HD_TL = Rates_HD(TL_filter);
        shifted_HD_TL = trueHD_TL - HD_PFDs(hd_cell);
        shifted_HD_TL(shifted_HD_TL>180) = shifted_HD_TL(shifted_HD_TL>180) - 360;
        hd_tuning_curve_TL = generate_tuning_curve(shifted_HD_TL,Rates_HD_TL,nBins_shift,binEdges_shift);
        hd_tuning_curve_TL(end) = [];
        if any(hd_tuning_curve_TL<0)
            disp('Something odd TL');
        end
        shifted_tuning_curves_TL(hd_cell,:) = hd_tuning_curve_TL;
        
          %ADN
        Rates_ADN_TL = Rates_ADN(TL_filter);
        adn_tuning_curve_TL = generate_tuning_curve(shifted_HD_TL,Rates_ADN_TL,nBins_shift,binEdges_shift);
        adn_tuning_curve_TL(end) = [];
        if any(adn_tuning_curve_TL<0)
            disp('Something odd TL');
        end
        shifted_tuning_curves_TL_ADN(hd_cell,:) = adn_tuning_curve_TL;
        
        %Bottom Left Quadrant
            %RSC
        Rates_HD_BL = Rates_HD(BL_filter);
        shifted_HD_BL = trueHD_BL - HD_PFDs(hd_cell);
        shifted_HD_BL(shifted_HD_BL>180) = shifted_HD_BL(shifted_HD_BL>180) - 360;
        hd_tuning_curve_BL = generate_tuning_curve(shifted_HD_BL,Rates_HD_BL,nBins_shift,binEdges_shift);
        hd_tuning_curve_BL(end) = [];
        if any(hd_tuning_curve_BL<0)
            disp('Something odd BL');
        end
        shifted_tuning_curves_BL(hd_cell,:) = hd_tuning_curve_BL;
        
          %ADN
        Rates_ADN_BL = Rates_ADN(BL_filter);
        adn_tuning_curve_BL = generate_tuning_curve(shifted_HD_BL,Rates_ADN_BL,nBins_shift,binEdges_shift);
        adn_tuning_curve_BL(end) = [];
        if any(adn_tuning_curve_BL<0)
            disp('Something odd BL');
        end
        shifted_tuning_curves_BL_ADN(hd_cell,:) = adn_tuning_curve_BL;
        
    end
    %Getting peaks and separations from autocorrelation
    autocorr_threshold = 0.5;
    
    [overall_autocorr_rschd(hd_cell,:), corrpks_all_rschd(hd_cell,1),corrpks_all_rschd(hd_cell,2:end)] = ...
        tuning_curve_autocorrelation(autocorr_threshold,hd_tuning_curve);
    
    if(sum(ADNRates(:,1)>0.01))
        [overall_autocorr_adn(hd_cell,:), corrpks_all_adn(hd_cell,1),corrpks_all_adn(hd_cell,2:end)] = ...
            tuning_curve_autocorrelation(autocorr_threshold,adn_tuning_curve);
    end
    
        if hd_cell == 24;
            plot_tuning_curve_v2(angles,compartments,0,'save',hd_cell,test_type,'HD',hd_tuning_curve);
            plot_tuning_curve_v2(angles,compartments,0,'save',hd_cell,test_type,'ADN',adn_tuning_curve);
        end
        
        if hd_cell == 124;
            plot_tuning_curve_v2(angles,compartments,0,'save',hd_cell,test_type,'HD',hd_tuning_curve);
            plot_tuning_curve_v2(angles,compartments,0,'save',hd_cell,test_type,'ADN',adn_tuning_curve);
        end
        
   
    
end


%% PLOT ALL THE FOUR QUADRANT STATISTICS
if strcmpi(test_type,'s')|| strcmpi(test_type,'d')

%calculate PVector of each tuning curve (in each quadrant) and then take stats from distribution

    shifts = -180:6:180;
    
    shifts_sin = sind(shifts); %sine
    shifts_cos = cosd(shifts); %cosine

    %RSC PVectors
    vector_1 = sum(bsxfun(@times,shifted_tuning_curves_TR,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,shifted_tuning_curves_TR,shifts_cos),2);
    TR_PV_all = atan2d(vector_1,vector_2);
    
    vector_1 = sum(bsxfun(@times,shifted_tuning_curves_BR,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,shifted_tuning_curves_BR,shifts_cos),2);
    BR_PV_all = atan2d(vector_1,vector_2);
    
    vector_1 = sum(bsxfun(@times,shifted_tuning_curves_TL,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,shifted_tuning_curves_TL,shifts_cos),2);
    TL_PV_all = atan2d(vector_1,vector_2);
    
    vector_1 = sum(bsxfun(@times,shifted_tuning_curves_BL,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,shifted_tuning_curves_BL,shifts_cos),2);
    BL_PV_all = atan2d(vector_1,vector_2);
    
    %ADN PVectors
    vector_1 = sum(bsxfun(@times,shifted_tuning_curves_TR_ADN,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,shifted_tuning_curves_TR_ADN,shifts_cos),2);
    TR_PV_all_ADN = atan2d(vector_1,vector_2);
    
    vector_1 = sum(bsxfun(@times,shifted_tuning_curves_BR_ADN,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,shifted_tuning_curves_BR_ADN,shifts_cos),2);
    BR_PV_all_ADN = atan2d(vector_1,vector_2);
    
    vector_1 = sum(bsxfun(@times,shifted_tuning_curves_TL_ADN,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,shifted_tuning_curves_TL_ADN,shifts_cos),2);
    TL_PV_all_ADN = atan2d(vector_1,vector_2);
    
    vector_1 = sum(bsxfun(@times,shifted_tuning_curves_BL_ADN,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,shifted_tuning_curves_BL_ADN,shifts_cos),2);
    BL_PV_all_ADN = atan2d(vector_1,vector_2);

    %Display means of the shifts
    TR_mean_all = mean(TR_PV_all)
    BR_mean_all = mean(BR_PV_all)
    TL_mean_all = mean(TL_PV_all)
    BL_mean_all = mean(BL_PV_all)

    TR_std_all = std(TR_PV_all)
    BR_std_all = std(BR_PV_all)
    TL_std_all = std(TL_PV_all)
    BL_std_all = std(BL_PV_all)
    
    TR_mean_all_ADN = mean(TR_PV_all_ADN)
    BR_mean_all_ADN = mean(BR_PV_all_ADN)
    TL_mean_all_ADN = mean(TL_PV_all_ADN)
    BL_mean_all_ADN = mean(BL_PV_all_ADN)

    TR_std_all_ADN = std(TR_PV_all_ADN)
    BR_std_all_ADN = std(BR_PV_all_ADN)
    TL_std_all_ADN = std(TL_PV_all_ADN)
    BL_std_all_ADN = std(BL_PV_all_ADN)


%now do the average tuning curve
%RSC
    mean_shifted_TR = mean(shifted_tuning_curves_TR,1);
    mean_shifted_BR = mean(shifted_tuning_curves_BR,1);
    mean_shifted_TL = mean(shifted_tuning_curves_TL,1);
    mean_shifted_BL = mean(shifted_tuning_curves_BL,1);
    
    std_shifted_TR = std(shifted_tuning_curves_TR,0,1);
    std_shifted_BR = std(shifted_tuning_curves_BR,0,1);
    std_shifted_TL = std(shifted_tuning_curves_TL,0,1);
    std_shifted_BL = std(shifted_tuning_curves_BL,0,1);
    
    shifts = -180:6:180;
    
    shifts_sin = sind(shifts); %sine
    shifts_cos = cosd(shifts); %cosine
    
    vector_1 = sum(bsxfun(@times,mean_shifted_TR,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,mean_shifted_TR,shifts_cos),2);
    TR_PV = atan2d(vector_1,vector_2);
    
    vector_1 = sum(bsxfun(@times,mean_shifted_BR,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,mean_shifted_BR,shifts_cos),2);
    BR_PV = atan2d(vector_1,vector_2);
    
    vector_1 = sum(bsxfun(@times,mean_shifted_TL,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,mean_shifted_TL,shifts_cos),2);
    TL_PV = atan2d(vector_1,vector_2);
    
    vector_1 = sum(bsxfun(@times,mean_shifted_BL,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,mean_shifted_BL,shifts_cos),2);
    BL_PV = atan2d(vector_1,vector_2);
    
    all = [mean_shifted_TR, mean_shifted_BR, mean_shifted_TL, mean_shifted_BL]; 
    figure()
    plot([shifts,shifts,shifts,shifts],all);
    y_limits = ylim;
    close(gcf);
    
    figure();
    subplot(2,2,1)
    [~,~] = boundedline(shifts,mean_shifted_TL,std_shifted_TL);
    title(['Top Left: ',num2str(TL_PV)]);
    ylabel('r [Hz]');
    xlabel('Shift (deg)');
    xlim([-180 180]);
    ylim(y_limits);
    set(gca,'Xtick',-180:60:180);
    hold on;
    line([TL_PV TL_PV],y_limits,'Linestyle','--');
    
    subplot(2,2,2);
    [~,~] = boundedline(shifts,mean_shifted_TR,std_shifted_TR);
    title(['Top Right: ',num2str(TR_PV)]);
    ylabel('r [Hz]');
    xlabel('Shift (deg)');
    xlim([-180 180]);
    ylim(y_limits);
    set(gca,'Xtick',-180:60:180);
    hold on;
    line([TR_PV TR_PV],y_limits,'Linestyle','--');
    
    subplot(2,2,3);
    [~,~] = boundedline(shifts,mean_shifted_BL,std_shifted_BL);
    title(['Bottom Left: ',num2str(BL_PV)]);
    ylabel('r [Hz]');
    xlabel('Shift (deg)');
    xlim([-180 180]);
    ylim(y_limits);
    set(gca,'Xtick',-180:60:180);
    hold on;
    line([BL_PV BL_PV],y_limits,'Linestyle','--');
    
    subplot(2,2,4);
    [~,~] = boundedline(shifts,mean_shifted_BR,std_shifted_BR);
    title(['Bottom Right: ',num2str(BR_PV)]);
    ylabel('r [Hz]');
    xlabel('Shift (deg)');
    xlim([-180 180]);
    ylim(y_limits);
    set(gca,'Xtick',-180:60:180);
    hold on;
    line([BR_PV BR_PV],y_limits,'Linestyle','--');
    
    if strcmpi(test_type,'s')
        saveas(gcf,'_shift_analysis_test','fig')
        close(gcf);
    else
        saveas(gcf,'_shift_analysis_disor','fig')
        close(gcf);
    end
    
    figure()
    plot(shifts,mean_shifted_TL,'Linewidth',3.0,'color',[0 0 1]);
    hold on
    plot(shifts,mean_shifted_TR,'Linewidth',3.0,'color',[1 0 0]);
    plot(shifts,mean_shifted_BL,'Linewidth',3.0,'color',[0 1 1]);
    plot(shifts,mean_shifted_BR,'Linewidth',3.0,'color',[1 0.5 0]);
    xlabel('Shift (deg)');
    xlim([-180 180]);
    set(gca,'Xtick',-180:60:180);
    
    if strcmpi(test_type,'s')
        saveas(gcf,'_shift_overlay_test','fig')
        close(gcf);
    else
        saveas(gcf,'_shift_overlay_disor','fig')
        close(gcf);
    end
    
    
%ADN
    mean_shifted_TR_ADN = mean(shifted_tuning_curves_TR_ADN,1);
    mean_shifted_BR_ADN = mean(shifted_tuning_curves_BR_ADN,1);
    mean_shifted_TL_ADN = mean(shifted_tuning_curves_TL_ADN,1);
    mean_shifted_BL_ADN = mean(shifted_tuning_curves_BL_ADN,1);
    
    std_shifted_TR_ADN = std(shifted_tuning_curves_TR_ADN,0,1);
    std_shifted_BR_ADN = std(shifted_tuning_curves_BR_ADN,0,1);
    std_shifted_TL_ADN = std(shifted_tuning_curves_TL_ADN,0,1);
    std_shifted_BL_ADN = std(shifted_tuning_curves_BL_ADN,0,1);
    
    shifts = -180:6:180;
    
    shifts_sin = sind(shifts); %sine
    shifts_cos = cosd(shifts); %cosine
    
    vector_1 = sum(bsxfun(@times,mean_shifted_TR_ADN,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,mean_shifted_TR_ADN,shifts_cos),2);
    TR_PV_ADN = atan2d(vector_1,vector_2);
    
    vector_1 = sum(bsxfun(@times,mean_shifted_BR_ADN,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,mean_shifted_BR_ADN,shifts_cos),2);
    BR_PV_ADN = atan2d(vector_1,vector_2);
    
    vector_1 = sum(bsxfun(@times,mean_shifted_TL_ADN,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,mean_shifted_TL_ADN,shifts_cos),2);
    TL_PV_ADN = atan2d(vector_1,vector_2);
    
    vector_1 = sum(bsxfun(@times,mean_shifted_BL_ADN,shifts_sin),2);
    vector_2 = sum(bsxfun(@times,mean_shifted_BL_ADN,shifts_cos),2);
    BL_PV_ADN = atan2d(vector_1,vector_2);
    
    all = [mean_shifted_TR_ADN, mean_shifted_BR_ADN, mean_shifted_TL_ADN, mean_shifted_BL_ADN]; 
    figure()
    plot([shifts,shifts,shifts,shifts],all);
    y_limits = ylim;
    close(gcf);
    
    figure();
    subplot(2,2,1)
    [~,~] = boundedline(shifts,mean_shifted_TL_ADN,std_shifted_TL_ADN);
    title(['Top Left: ',num2str(TL_PV_ADN)]);
    ylabel('r [Hz]');
    xlabel('Shift (deg)');
    xlim([-180 180]);
    ylim(y_limits);
    set(gca,'Xtick',-180:60:180);
    hold on;
    line([TL_PV_ADN TL_PV_ADN],y_limits,'Linestyle','--');
    
    subplot(2,2,2);
    [~,~] = boundedline(shifts,mean_shifted_TR_ADN,std_shifted_TR_ADN);
    title(['Top Right: ',num2str(TR_PV_ADN)]);
    ylabel('r [Hz]');
    xlabel('Shift (deg)');
    xlim([-180 180]);
    ylim(y_limits);
    set(gca,'Xtick',-180:60:180);
    hold on;
    line([TR_PV_ADN TR_PV_ADN],y_limits,'Linestyle','--');
    
    subplot(2,2,3);
    [~,~] = boundedline(shifts,mean_shifted_BL_ADN,std_shifted_BL_ADN);
    title(['Bottom Left: ',num2str(BL_PV_ADN)]);
    ylabel('r [Hz]');
    xlabel('Shift (deg)');
    xlim([-180 180]);
    ylim(y_limits);
    set(gca,'Xtick',-180:60:180);
    hold on;
    line([BL_PV_ADN BL_PV_ADN],y_limits,'Linestyle','--');
    
    subplot(2,2,4);
    [~,~] = boundedline(shifts,mean_shifted_BR_ADN,std_shifted_BR_ADN);
    title(['Bottom Right: ',num2str(BR_PV_ADN)]);
    ylabel('r [Hz]');
    xlabel('Shift (deg)');
    xlim([-180 180]);
    ylim(y_limits);
    set(gca,'Xtick',-180:60:180);
    hold on;
    line([BR_PV_ADN BR_PV_ADN],y_limits,'Linestyle','--');
    
     if strcmpi(test_type,'s')
    saveas(gcf,'_shift_analysis_test_ADN','fig')
    close(gcf);
     else 
          saveas(gcf,'_shift_analysis_disor_ADN','fig')
    close(gcf);
     end
    
    figure()
    plot(shifts,mean_shifted_TL_ADN,'Linewidth',3.0,'color',[0 0 1]);
    hold on
    plot(shifts,mean_shifted_TR_ADN,'Linewidth',3.0,'color',[1 0 0]);
    plot(shifts,mean_shifted_BL_ADN,'Linewidth',3.0,'color',[0 1 1]);
    plot(shifts,mean_shifted_BR_ADN,'Linewidth',3.0,'color',[1 0.5 0]);
    xlabel('Shift (deg)');
    xlim([-180 180]);
    set(gca,'Xtick',-180:60:180);
    
     
     if strcmpi(test_type,'s')
    saveas(gcf,'_shift_overlay_test_ADN','fig')
    close(gcf);
     else 
          saveas(gcf,'_shift_overlay_disor_ADN','fig')
    close(gcf);
     end
end

%% NOW PLOT THE RSC HD SUMMARIES

%Look at average autocorrelations for all HD cells
figure()
[~,~] = boundedline(angles(1:60),mean(overall_autocorr_rschd(:,1:60)),std(overall_autocorr_rschd(:,1:60)));
xlim([0 2*pi]);
set(gca,'XTick',0:(2/6)*pi:2*pi);
set(gca,'XTickLabel',0:60:360);
ylim([0 3]);
set(gca,'Fontsize',14);
title('All HD Average');
xlabel('Rotation');
ylabel('Correlation');
if (strcmpi(save_mode, 'save'))
    if strcmpi(test_type,'s')
        saveas(gcf,'_RSCHD_autocorr_summary_test', 'fig');
    elseif strcmpi(test_type,'d')
        saveas(gcf,'_RSCHD_autocorr_summary_disor', 'fig');
    else
        saveas(gcf,'_RSCHD_autocorr_summary', 'fig');
    end
    close(gcf);
    
end

%Look at autocorrelation heat map for all HD cells
figure()
cmax = max(max(overall_autocorr_rschd));
imagesc(overall_autocorr_rschd);
set(gca,'YDir','normal')
set(gca,'YTickLabel',hd_cells/10:hd_cells/10:hd_cells);
set(gca,'XTick',1:10:61);
set(gca,'XTickLabel',0:60:360);
set(gca,'TickLength',[ 0 0 ]);
caxis([0 cmax]);
set(gca,'Fontsize',14);
title('All HD');
if (strcmpi(save_mode, 'save'))
    if strcmpi(test_type,'s')
        saveas(gcf,'_RSCHD_autocorr_test', 'fig');
    elseif strcmpi(test_type,'d')
        saveas(gcf,'_RSCHD_autocorr_disor', 'fig');
    else
        saveas(gcf,'_RSCHD_autocorr', 'fig');
    end
    close(gcf);
    
end

%% NOW PLOT THE ADN HD SUMMARIES

%Look at average autocorrelations for all ADN HD cells
figure()
[~,~] = boundedline(angles(1:60),mean(overall_autocorr_adn(:,1:60)),std(overall_autocorr_adn(:,1:60)));
xlim([0 2*pi]);
set(gca,'XTick',0:(2/6)*pi:2*pi);
set(gca,'XTickLabel',0:60:360);
ylim([0 14])
set(gca,'Fontsize',14);
title('All HD Average');
xlabel('Rotation');
ylabel('Correlation');
if (strcmpi(save_mode, 'save'))
    if strcmpi(test_type,'s')
        saveas(gcf,'_ADNHD_autocorr_summary_test', 'fig');
    elseif strcmpi(test_type,'d')
        saveas(gcf,'_ADNHD_autocorr_summary_disor', 'fig');
    else
        saveas(gcf,'_ADNHD_autocorr_summary', 'fig');
    end
    close(gcf);
    
end
%Look at autocorrelation heat map for all ADN HD cells
figure()
cmax = max(max(overall_autocorr_adn));
imagesc(overall_autocorr_adn);
set(gca,'YDir','normal')
set(gca,'YTickLabel',hd_cells/10:hd_cells/10:hd_cells);
set(gca,'XTick',1:10:61);
set(gca,'XTickLabel',0:60:360);
set(gca,'TickLength',[ 0 0 ]);
caxis([0 cmax]);
set(gca,'Fontsize',14);
title('All HD');
if (strcmpi(save_mode, 'save'))
    if strcmpi(test_type,'s')
        saveas(gcf,'_ADNHD_autocorr_test', 'fig');
    elseif strcmpi(test_type,'d')
        saveas(gcf,'_ADNHD_autocorr_disor', 'fig');
    else
        saveas(gcf,'_ADNHD_autocorr', 'fig');
    end
    close(gcf);
    
end

%% NOW SAVE THE DATA

if strcmpi(test_type,'s')
    %RSC HD
    save('overall_tuning_curve_rschd_test.mat','overall_tuning_curve_rschd');
    save('overall_autocorr_rschd_test.mat','overall_autocorr_rschd');
    save('corrpks_all_rschd_test.mat','corrpks_all_rschd');
    
    %ADN HD
    save('overall_tuning_curve_adn_test.mat','overall_tuning_curve_adn');
    save('overall_autocorr_adn_test.mat','overall_autocorr_adn');
    save('corrpks_all_adn_test.mat','corrpks_all_adn');
elseif strcmpi(test_type,'d')
    %RSC HD
    save('overall_tuning_curve_rschd_disor.mat','overall_tuning_curve_rschd');
    save('overall_autocorr_rschd_disor.mat','overall_autocorr_rschd');
    save('corrpks_all_rschd_disor.mat','corrpks_all_rschd');
    
    %ADN HD
    save('overall_tuning_curve_adn_disor.mat','overall_tuning_curve_adn');
    save('overall_autocorr_adn_disor.mat','overall_autocorr_adn');
    save('corrpks_all_adn_disor.mat','corrpks_all_adn');
else
    %RSC HD
    save('overall_tuning_curve_rschd.mat','overall_tuning_curve_rschd');
    save('overall_autocorr_rschd.mat','overall_autocorr_rschd');
    save('corrpks_all_rschd.mat','corrpks_all_rschd');
    
    %ADN HD
    save('overall_tuning_curve_adn.mat','overall_tuning_curve_adn');
    save('overall_autocorr_adn.mat','overall_autocorr_adn');
    save('corrpks_all_adn.mat','corrpks_all_adn');
end

end
function rates = file_load(cells, steps, fname)

rates = zeros(cells, steps);

fid = fopen(fname, 'rb');

rates = fread(fid, [steps, cells], 'float32')';

fclose(fid);


end

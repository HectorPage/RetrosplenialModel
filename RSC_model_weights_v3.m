function RSC_model_weights_v2(mode,epochs,varargin)
%plots weights for a RSC/conjunctive cell model

%Created by and copyright held by Dr. Hector Page 23/09/16

%USAGE:
% Make the first argument either 'file' to load in weight profiles from
%file, or 'simulation' to pass weights directly from simulation

%If mode = 'file', make useage fullySOWeights('file', epochs, varargin), with varagins as follows:
                                                % conjunctive_cells = vargin{1};
                                                % hd_cells  = vargin{2};    
                                                % vis_cells = vargin{3};
                                                % adn_cells - vargin{3};
                                        

%If mode = 'simluation', make usage as follows:
                                                % VIStoCONJ_weights_init  = vargin{1};
                                                % VIStoCONJ_weights_fin   = vargin{2};
    
                                                % HDtoCONJ_weights_init   = vargin{3};
                                                % HDtoCONJ_weights_fin    = vargin{4};

                                                % CONJtoHD_weights_init   = vargin{5};
                                                % CONJtoHD_weights_fin    = vargin{6};
                                                 
                                                % VIStoCONJ_weight_change = vargin{7};
                                                % HDtoCONJ_weight_change  = vargin{8};
                                                % CONJtoHD_weight_change  = vargin{9};
                                                
                                                %ADNtoHD_weights = varargin{10};
                                                %ADNtoADN_weights = vargin{12}
                                                %HDtoADN_weights = vargin{12};
                                                
                                                

if strcmpi(mode,'file')

    conjunctive_cells = varargin{1}; %make sure vargins are correctly ordered....
    hd_cells = varargin{2};
    vis_cells = varargin{3};
    adn_cells = varargin{4};
    
    % LOADING IN WEIGHTS - CHECK DIMENSIONS AND FORMAT ARE CORRECT!
    fid = fopen('vis_to_conjunctive_weights_initial.bdat', 'rb');
    VIStoCONJ_weights_init = fread(fid, [vis_cells,conjunctive_cells], 'float32')';
    fclose(fid);

    fid = fopen('vis_to_conjunctive_weights_final.bdat', 'rb');
    VIStoCONJ_weights_fin = fread(fid, [vis_cells,conjunctive_cells], 'float32')';
    fclose(fid);

    fid = fopen('HD_to_conjunctive_weights_initial.bdat', 'rb');
    HDtoCONJ_weights_init = fread(fid, [hd_cells,conjunctive_cells], 'float32')';
    fclose(fid);

    fid = fopen('conjunctive_to_HD_weights_initial.bdat', 'rb');
    CONJtoHD_weights_init = fread(fid, [conjunctive_cells,hd_cells], 'float32')';
    fclose(fid);

    fid = fopen('HD_to_conjunctive_weights_final.bdat', 'rb');
    HDtoCONJ_weights_fin = fread(fid, [hd_cells,conjunctive_cells], 'float32')';
    fclose(fid);

    fid = fopen('conjunctive_to_HD_weights_final.bdat', 'rb');
    CONJtoHD_weights_fin = fread(fid, [conjunctive_cells,hd_cells], 'float32')';
    fclose(fid);
    
    fid = fopen('ADN_to_HD_weights.bdat', 'rb');
    ADNtoHD_weights = fread(fid, [adn_cells,hd_cells], 'float32')';
    fclose(fid);
    
    fid = fopen('ADN_to_ADN_weights.bdat', 'rb');
    ADNtoADN_weights = fread(fid, [adn_cells,adn_cells], 'float32')';
    fclose(fid);
    
    fid = fopen('HD_to_ADN_weights.bdat', 'rb');
    HDtoADN_weights = fread(fid, [hd_cells,adn_cells], 'float32')';
    fclose(fid);

    %LOADING IN WEIGHT CHANGES
    fileID = fopen('VIStoCONJ_weight_change.bdat','rb');
    VIStoCONJ_weight_change = fread(fid, [1,epochs/100], 'float32')';
    fclose(fileID);
    
    fileID = fopen('HDtoCONJ_weight_change.bdat','rb');
    HDtoCONJ_weight_change  = fread(fid, [1,epochs/100], 'float32')';
    fclose(fileID);
    
    fileID = fopen('CONJtoHD_weight_change.bdat','rb');
    CONJtoHD_weight_change  = fread(fid, [1,epochs/100], 'float32')';
    fclose(fileID);


elseif strcmpi(mode, 'simulation')
    
    VIStoCONJ_weights_init  = varargin{1};
    VIStoCONJ_weights_fin   = varargin{2};
    
    HDtoCONJ_weights_init   = varargin{3};
    HDtoCONJ_weights_fin    = varargin{4};
    
    CONJtoHD_weights_init   = varargin{5};
    CONJtoHD_weights_fin    = varargin{6};
    
    VIStoCONJ_weight_change = varargin{7};
    HDtoCONJ_weight_change  = varargin{8};
    CONJtoHD_weight_change  = varargin{9};   
    
    ADNtoHD_weights = varargin{10};
    ADNtoADN_weights = varargin{11};
    HDtoADN_weights = varargin{12};
    
else
    error('First argument must be ''file'' or ''simulation''!');
end

%% Weight changes
figure()
subplot(3,1,1);
plot(VIStoCONJ_weight_change,'Linewidth',2.0);
title('Vis to Conj Weight Change');
subplot(3,1,2);
plot(HDtoCONJ_weight_change,'Linewidth',2.0);
title('HD to Conj Weight Change');
subplot(3,1,3);
plot(CONJtoHD_weight_change,'Linewidth',2.0);
title('Conj to HD Weight Change');
savefig('_WeightChanges');
close(gcf);

%% VIS to conjunctive weights

figure()
subplot(2,1,1)
map = imagesc(VIStoCONJ_weights_init);
cbar = [0 max(max(VIStoCONJ_weights_init))];
caxis(cbar);
colormap(jet);
cb = colorbar;
ylabel('Vis Cell');
xlabel('Conjunctive Cell');
ylabel(cb, 'Weight', 'Fontsize',12)
title('Initial Vis to Conjunctive Weights')

subplot(2,1,2)
map = imagesc(VIStoCONJ_weights_fin);
cbar = [0 max(max(VIStoCONJ_weights_fin))];
caxis(cbar);
colormap(jet);
cb = colorbar;
ylabel('Vis Cell');
xlabel('Conjunctive Cell');
ylabel(cb, 'Weight', 'Fontsize',12)
title('Final Vis to Conjunctive Weights')
savefig('_VistoCONJWeights');
close(gcf);

%% INITIAL HD to conjunctive weights
figure()
subplot(2,1,1)
map = imagesc(HDtoCONJ_weights_init);
cbar = [0 max(max(HDtoCONJ_weights_init))];
caxis(cbar);
colormap(jet);
cb = colorbar;
ylabel('RSC HD Cell');
xlabel('Conjunctive Cell');
ylabel(cb, 'Weight', 'Fontsize',12);
title('Initial HD to Conjunctive Weights');

subplot(2,1,2)

%% FINAL HD to conjunctive weights
map = imagesc(HDtoCONJ_weights_fin);
cbar = [0 max(max(HDtoCONJ_weights_fin))];
caxis(cbar);
colormap(jet);
cb = colorbar;
ylabel('RSC HD Cell');
xlabel('Conjunctive Cell');
ylabel(cb, 'Weight', 'Fontsize',12);
title('Final HD to conjunctive Weights');
savefig('_HDtoCONJWeights');
close(gcf);


%% INITIAL conjunctive to HD weights
figure()
subplot(2,1,1)
map = imagesc(CONJtoHD_weights_init);
cbar = [0 max(max(CONJtoHD_weights_init))];
caxis(cbar);
colormap(jet);
cb = colorbar;
ylabel('Conjunctive Cell');
xlabel('RSC HD Cell');
ylabel(cb, 'Weight', 'Fontsize',12);
title('Initial conjunctive to HD Weights');

%% FINAL conjunctive to HD weights

subplot(2,1,2)
map = imagesc(CONJtoHD_weights_fin);
cbar = [0 max(max(CONJtoHD_weights_fin))];
caxis(cbar);
colormap(jet);
cb = colorbar;
ylabel('Conjunctive Cell');
xlabel('RSC HD Cell');
ylabel(cb, 'Weight', 'Fontsize',12);
title('Final conjunctive to HD Weights');
savefig('_CONJtoHDWeights');
close(gcf);

%% Weights to and from ADN
figure()
subplot(3,1,1);
map = imagesc(ADNtoHD_weights);
cbar = [0 max(max(ADNtoHD_weights))];
caxis(cbar);
colormap(jet);
cb = colorbar;
ylabel('HD Cell');
xlabel('ADN Cell');
subplot(3,1,2);
map = imagesc(HDtoADN_weights);
cbar = [0 max(max(HDtoADN_weights))];
caxis(cbar);
colormap(jet);
cb = colorbar;
ylabel('ADN Cell');
xlabel('HD Cell');
subplot(3,1,3);
map = imagesc(ADNtoADN_weights);
cbar = [0 max(max(ADNtoADN_weights))];
caxis(cbar);
colormap(jet);
cb = colorbar;
ylabel('ADN Cell');
xlabel('ADN Cell');
savefig('_ADNtofromallWeights');
close(gcf);


end



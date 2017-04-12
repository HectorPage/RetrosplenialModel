function [time, training_time, HD_cells,ADN_cells,conjunctive_cells, num_landmarks, phi_HDconjunctive, phi_conjunctiveHD, phi_vis,...
phi_ADNHD,phi_HDADN,phi_ADNADN,tau_HD, tau_conjunctive, tau_ADN, timestep_size, HD_inhibition, conjunctive_inhibition, ADN_inhibition,...
conduction_delay,learning_rate, vis_conjunctive_weight_sigma, HD_conjunctive_weight_sigma, conjunctive_HD_weight_sigma, vis_sigma,...
ADN_ADN_weight_sigma, ADN_HD_weight_sigma, HD_ADN_weight_sigma,PI_noise_percentage,PI_noise_width, PI_strength, PI_sigma, LTD] = read_RSC_v2_params(filename)
%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/hector.page/flip_cell_model/_MATLAB_MODEL/_bipolar_params.xlsx
%    Worksheet: Sheet1
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2016/10/13 10:30:30

%% Import the data
[~, ~, raw] = xlsread(filename,'Sheet1');
raw = raw(:,end);

%% Create output variable
params = reshape([raw{:}],size(raw));

time = params(1);
training_time  = params(2);
HD_cells  = params(3);
ADN_cells  = params(4);
conjunctive_cells  = params(5);
num_landmarks  = params(6);
phi_HDconjunctive  = params(7);
phi_conjunctiveHD  = params(8);
phi_vis  = params(9);
phi_ADNHD  = params(10);
phi_HDADN  = params(11);
phi_ADNADN  = params(12);
tau_HD  = params(13);
tau_conjunctive  = params(14);
tau_ADN  = params(15);
timestep_size  = params(16);
HD_inhibition  = params(17);
conjunctive_inhibition  = params(18);
ADN_inhibition  = params(19);
conduction_delay  = params(20);
learning_rate  = params(21);
vis_conjunctive_weight_sigma  = params(22);
HD_conjunctive_weight_sigma  = params(23);
conjunctive_HD_weight_sigma  = params(24);
vis_sigma  = params(25);
ADN_ADN_weight_sigma = params(26);
ADN_HD_weight_sigma = params(27);
HD_ADN_weight_sigma = params(28);
PI_noise_percentage = params(29);
PI_noise_width = params(30);
PI_strength = params(31);
PI_sigma = params(32);
LTD = params(33);
 

%% Clear temporary variables
clearvars raw params;

end


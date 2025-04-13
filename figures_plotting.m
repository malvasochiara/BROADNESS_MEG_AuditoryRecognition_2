%%
%% FIGURE 2A: 2 FIGURES, ONE FOR EACH COMPONENT, SHOWING ALL THE CONDITIONS IN EACH FIGURE
% loading data
clear all
addpath('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure2\result02_fdr.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure2\timeserie_PCAonmaineffect.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time

col_l = 1 ; %1 for significant time-windows with different colors; 0 = only grey color
ylimm = []; %amplitude limits; leave empty [] for automatic adjustment

PCs = 2;  % Number of principal components you want to plot

export_l = 0; %1 = export images; 0 = not
% cnt = 0;
data = permute(J,[2 1 4 3]); 
close all
for cc = 1:PCs %over principal components
    lineplotdum = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
    lineplot = lineplotdum; 
    S = [];
    S.sbrim = 0.15;
    S.ii = cc;
    S.conds = {'Old ','NewT1','NewT2','NewT3','NewT4'};
    S.data = data(cc,:,:,:);
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.time_real = time_sel(1:1025);
    S.colorline = [1 0 0; 0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];        
    if export_l == 1
        S.legendl = 0;
    else
        S.legendl = 1;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = 1;
    S.condition_n = 1:5; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
    S.ROIs_labels = {1, 2}; 
    S.lineplot = lineplot;
    bum = [];
    %%%SIGNIFICATIVITA%%%
    signtp_col = [];
    for ss = 1:4 %over contrasts between conditions
        sbam = size(fdr_results{cc,ss},2);
        clear signtp_temp;
        for ll = 1:sbam %over the number of the significant time-windows for contrast ss
            bum = cat(2,bum,{fdr_results{cc,ss}(ll).time_interval}); %concatenating all significant time-windows, for one contrast at a time
            
            signtp_temp(ll) = (ss+1); %getting the color code (number of the significant time-windows for each contrast)
        end
        signtp_col = cat(2, signtp_col, signtp_temp);
    end
    S.signtp = bum;
    if col_l == 1
        S.signtp_col = signtp_col;
        S.colorsign = S.colorline;
    else
        S.signtp_col = [];
        S.colorsign = {'*';'+'};
    end
    
    S.groups = {1}; % this is one bc you subjects are not divided in groups, set this to 1 and don't ever touch it 
    S.gsubj = {1:83}; %same as above, this contains the indices of subjects and they range from 0 to 83 since your subjects are not divided in groups
    waveplot_groups_local_v2(S) %actual function
    
    if export_l == 1
        exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure2\PCAmaineffect_component_0' num2str(cc) '.pdf'],'Resolution',300)
    end
    %end
    hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 2B: 5 FIGURES, ONE FOR EACH CONDITION, SHOWING ALL THE COMPONENTS IN EACH FIGURE
clear all;
% loading data
addpath('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure2\result02_fdr.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure2\timeserie_PCAonmaineffect.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time

col_l = 0; %1 for significant time-windows with different colors; 1 = only grey color
ylimm = [-1900 800]; %amplitude limits; leave empty [] for automatic adjustment

conds = 5;  % number of conditions

export_l = 1; %1 = export images; 0 = not

data = permute(J,[3 1 4 2]);
close all
for cc = 1:conds %over conditions
    lineplotdum = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
    lineplot = repmat(lineplotdum,cc,1);
    S = [];
    S.sbrim = 0.15;
    S.ii = cc;

    S.conds = {'PC 1 ','PC 2'};
    S.data = data(cc,:,:,:);
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.time_real = time_sel(1:1025);
    S.colorline = [0.7 0.2 0.2;0.15 0.25 0.6];        
    if export_l == 1
        S.legendl = 0;
    else
        S.legendl = 1;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = 1;
    S.condition_n = 1:2; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
    S.ROIs_labels = {1, 2, 3, 4, 5}; 
    S.lineplot = lineplot;
    S.signtp = [];
    S.signtp_col = [];
    S.groups = {1}; % this is one bc you subjects are not divided in groups, set this to 1 and don't ever touch it 
    S.gsubj = {1:83}; %same as above, this containes the indices os subjects and they range from 0 to 83 since your subjects are not divided in groups
    waveplot_groups_local_v2(S) %actual function
    if export_l == 1
        exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure2\PCAmaineffect_condition_0' num2str(cc) '.pdf'],'Resolution',300)
    end
    hold on;
end

%% FIGURE 3 VAR EXP
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VAR EXPLAINED  FOR THE TWO RANDOMIZATIONS
var_rand1 = load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure3_newcode\randomiazionoverTIME_varexp.mat');
var_rand2 = load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure3_newcode\randomiazionoverSPACELABELS_varexp.mat');
var_ref = load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure3_newcode\varexp_from_maineffect.mat'); % (?)

export_l = 1; % 1 to export the fig
colorline = [1 0 0; 0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];
figure;
plot(var_ref.vare(1:20), 'Color',colorline(1,:),'Linewidth', 3, 'DisplayName', 'No rand');
hold on;
plot(var_rand1.varef_averaged(1:20), 'Color',colorline(2,:),'Linewidth', 1.5, 'DisplayName', 'Rand over time');
hold on;
plot(var_rand2.varef_averaged(1:20), 'Color',colorline(3,:), 'Linewidth', 1.5 ,'DisplayName', 'Rand over space labels');
hold on;
grid minor
set(gcf,'color','w')
box on

if export_l == 0
    legend('show');
end

if export_l ==1
    exportgraphics(gcf,'C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure3_newcode\Randomization_varexplained.pdf','Resolution',300)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 3 TIME SERIE RAND OVER TIME RECONSTRUCTED ON RAND MATRIX
clear all
%loading data:
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time
addpath('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure3_newcode\randomiazionoverTIME_result02_fdr.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure3_newcode\randomiazionoverTIME_tsrand_averaged.mat');

col_l = 1 ; %1 for significant time-windows with different colors; 0 = only grey color
ylimm = []; %amplitude limits; leave empty [] for automatic adjustment

PCs = 2;  %IN YOUR CASE THESE ARE THE COMPONENTS, NOT THE REGION OF INTEREST

export_l = 1; %1 = export images; 0 = not

data = permute(J_rand_averaged,[2 1 3 4]);

% demeaning
data_demeaned = zeros(size(data,1),size(data,2),size(data,3),size(data,4));
for cc = 1:size(data,1) % over components
    for sub = 1:size(data,3) % over subjects
        for cond = 1: size(data,4) %over conditions
            meandata = mean(data(cc,:,sub,cond),2);
            data_demeaned(cc,:,sub,cond) = data(cc,:,sub,cond) - meandata;
        end
    end
end
close all
for cc = 1:PCs %over principal components
    lineplotdum = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
    lineplot = lineplotdum; 
    S = [];
    S.sbrim = 0.15;
    S.ii = cc;
    S.conds = {'Old ','NewT1','NewT2','NewT3','NewT4'};
    S.data = data_demeaned(cc,:,:,:);
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.time_real = time_sel(1:1025);
    S.colorline = [1 0 0; 0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];        
    if export_l == 1
    S.legendl = 0;
    else
    S.legendl = 1;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = 1;
    S.condition_n = 1:5; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
    S.ROIs_labels = {1, 2}; 
    S.lineplot = lineplot;
    bum = [];
    %%%SIGNIFICATIVITA%%%
    signtp_col = [];
    for ss = 1:4 %over contrasts between conditions
    sbam = size(fdr_results{cc,ss},2);
    if sbam == 0
        continue;
    end
    clear signtp_temp;
    for ll = 1:sbam %over the number of the significant time-windows for contrast ss
        bum = cat(2,bum,{fdr_results{cc,ss}(ll).time_interval}); %concatenating all significant time-windows, for one contrast at a time
        
        signtp_temp(ll) = (ss+1); %getting the color code (number of the significant time-windows for each contrast)
    end
    signtp_col = cat(2, signtp_col, signtp_temp);
    end
    S.signtp = bum;
    if col_l == 1
        S.signtp_col = signtp_col;
        S.colorsign = S.colorline;
    else
    S.signtp_col = [];
    S.colorsign = [];
    end
    
    S.groups = {1}; % this is one bc you subjects are not divided in groups, set this to 1 and don't ever touch it 
    S.gsubj = {1:83}; %same as above, this containes the indices os subjects and they range from 0 to 83 since your subjects are not divided in groups
    waveplot_groups_local_v2(S) %actual function
    %title(lab(ROIs_to_AAL{ii,1}(pp),:))
    
    if export_l == 1
        exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure3_newcode\RandomizationTIME_timeserie_component_0' num2str(cc) '.pdf'],'Resolution',300)
    end
    hold on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 3 TIME SERIE RAND OVER SPACE LABELS RECONSTRUCTED ON RAND MATRIX

clear all
%loading data:
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time
addpath('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure3_newcode\randomiazionoverSPACELABELS_result02_fdr.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure3_newcode\randomiazionoverSPACELABELS_tsrand.mat');
% J = zeros(1025,2,1,5);
% J(:,:,1,:) = J_rand(:,:,:);
col_l = 1 ; %1 for significant time-windows with different colors; 0 = only grey color
ylimm = []; %amplitude limits; leave empty [] for automatic adjustment

PCs = 2;  %IN YOUR CASE THESE ARE THE COMPONENTS, NOT THE REGION OF INTEREST

export_l = 1; %1 = export images; 0 = not

data =permute(J_rand,[2 1 3 4]);
close all
for cc = 1:PCs %over principal components
    lineplotdum = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
    lineplot = lineplotdum; 
    S = [];
    S.sbrim = 0.15;
    S.ii = cc;
    S.conds = {'Old ','NewT1','NewT2','NewT3','NewT4'};
    S.data = data(cc,:,:,:);
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.time_real = time_sel(1:1025);
    S.colorline = [1 0 0; 0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];        
    if export_l == 1
    S.legendl = 0;
    else
    S.legendl = 1;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = 1;
    S.condition_n = 1:5; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
    S.ROIs_labels = {1, 2}; 
    S.lineplot = lineplot;
    bum = [];
    %%%SIGNIFICATIVITA%%%
    signtp_col = [];
    for ss = 1:4 %over contrasts between conditions
    sbam = size(fdr_results{cc,ss},2);
    if sbam == 0
        continue;
    end
    clear signtp_temp;
    for ll = 1:sbam %over the number of the significant time-windows for contrast ss
        bum = cat(2,bum,{fdr_results{cc,ss}(ll).time_interval}); %concatenating all significant time-windows, for one contrast at a time
        
        signtp_temp(ll) = (ss+1); %getting the color code (number of the significant time-windows for each contrast)
    end
    signtp_col = cat(2, signtp_col, signtp_temp);
    end
    S.signtp = bum;
    if col_l == 1
        S.signtp_col = signtp_col;
        S.colorsign = S.colorline;
    else
    S.signtp_col = [];
    S.colorsign = [];
    end
    
    S.groups = {1}; % this is one bc you subjects are not divided in groups, set this to 1 and don't ever touch it 
    S.gsubj = {1:83}; %same as above, this containes the indices os subjects and they range from 0 to 83 since your subjects are not divided in groups
    waveplot_groups_local_v2(S) %actual function
    %title(lab(ROIs_to_AAL{ii,1}(pp),:))
    
    if export_l == 1
        exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\Figure3_newcode\RandomizationSPACELABELS_timeserie_component_0' num2str(cc) '.pdf'],'Resolution',300)
    end
    hold on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FIGURE S1: 6 FIGURES, ONE FOR EACH COMPONENT, SHOWING ALL THE CONDITIONS IN EACH FIGURE
clear all
% loading data
addpath('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS1\result02_fdr.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS1\timeserie_wcoeff_from_maineffect.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time

col_l = 1 ; %1 for significant time-windows with different colors; 0 = only grey color
ylimm = [-1950, 800]; %amplitude limits; leave empty [] for automatic adjustment

PCs = 6;  % number of principal components you want to plot

export_l = 0; %1 = export images; 0 = not
data = permute(J,[2 1 4 3]);
close all
for cc = 1:PCs %over original parcels (+1 which is the voxels which did not belong to any parcel)
    lineplotdum = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
    lineplot = lineplotdum; 
    S = [];
    S.sbrim = 0.15;
    S.ii = cc;
    S.conds = {'Old ','NewT1','NewT2','NewT3','NewT4'};
    S.data = data(cc,:,:,:);
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.time_real = time_sel(1:1025);
    S.colorline = [1 0 0; 0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];        
    if export_l == 1
        S.legendl = 0;
    else
        S.legendl = 1;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = 1;
    S.condition_n = 1:5; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
    S.ROIs_labels = {1, 2, 3, 4, 5, 6}; 
    S.lineplot = lineplot;
    bum = [];
    %%%SIGNIFICATIVITA%%%
    signtp_col = [];
    for ss = 1:4 %over contrasts between conditions
        sbam = size(fdr_results{cc,ss},2);
        if sbam == 0
            continue;
        end
        clear signtp_temp;
        for ll = 1:sbam %over the number of the significant time-windows for contrast ss
            bum = cat(2,bum,{fdr_results{cc,ss}(ll).time_interval}); %concatenating all significant time-windows, for one contrast at a time
            
            signtp_temp(ll) = (ss+1); %getting the color code (number of the significant time-windows for each contrast)
        end
        signtp_col = cat(2, signtp_col, signtp_temp);
    end
    S.signtp = bum;
    if col_l == 1
        S.signtp_col = signtp_col;
        S.colorsign = S.colorline;
    else
        S.signtp_col = [];
        S.colorsign = {'*';'+'};
    end
    
    S.groups = {1}; % this is one bc you subjects are not divided in groups, set this to 1 and don't ever touch it 
    S.gsubj = {1:83}; %same as above, this containes the indices os subjects and they range from 0 to 83 since your subjects are not divided in groups
    waveplot_groups_local_v2(S) %actual function
    
    if export_l == 1
        exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS1\PCAmaineffect_component_0' num2str(cc) '.pdf'],'Resolution',300)
    end
    %end
    hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE S2: EXPLAINED VARIANCE
clear all
%%%% EXPLAINED VAR %%%%
var_singlecond = load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\figureS2\singleconditions_varexp_from_maineffect.mat');
var_avcond = load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\figureS2\averagedconditions_varexp_from_maineffect.mat');


export_l = 0; % 1 to export the figure
colorline = [0.9 0.4 0.3;0.5 0.4 0.75;0.3 0.6 0.9;0.12 0.18 0.6;0.1 0.1 0.1; 0.7 0.3 0.5];
figure;
plot(var_avcond.varexp(1:20), 'Color',colorline(1,:),'Linewidth', 1.5, 'DisplayName', 'Averaged conditions');
hold on;
for cc = 1:size(var_singlecond.varexp,2)
    plot(var_singlecond.varexp(1:20,cc), 'Color',colorline(cc+1,:),'Linewidth', 1.5, 'DisplayName', ['Condition ' num2str(cc)]);
    hold on;
end
grid minor
set(gcf,'color','w')
box on

if export_l == 0
    legend('show');
end

if export_l ==1
    exportgraphics(gcf,'C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS2\singleVSaveragedcondition_varexplained.pdf','Resolution',300)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE S2: 5 FIGURES, ONE FOR EACH CONDITION, SHOWING ALL THE COMPONENTS IN EACH FIGURE FOR AVERAGED COND AND SINGLE COND 
clear all
ts_avcond = load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\figureS2\averagedconditions_timeserie_wcoeff_from_maineffect.mat');
ts_singlecond = load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\figureS2\singleconditions_timeserie_wcoeff_from_maineffect.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time

data = permute(cat(2,ts_avcond.J,ts_singlecond.J), [3 1 4 2]);

col_l = 0; %1 for significant time-windows with different colors; 1 = only grey color
ylimm = []; %amplitude limits; leave empty [] for automatic adjustment

conds = 5;  % IN YOUR CASE THESE ARE THE CONDITIONS, NOT 

export_l = 0; %1 = export images; 0 = not
close all
for cc = 1:conds %over original parcels (+1 which is the voxels which did not belong to any parcel)
    lineplotdum = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
    lineplot = repmat(lineplotdum,cc,1);
    S = [];
    S.sbrim = 0.15;
    S.ii = cc;

    S.conds = {'Averaged cond, PC 1','Averaged cond, PC 2', 'Single cond, PC 1', 'Single cond, PC 2'};

    S.data = data(cc,:,:,:);
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.time_real = time_sel(1:1025);
    S.colorline = [0.58 0 0;0.0784 0.1569 0.5255;1.0 0.01 0.0;0.1882 0.4902 0.8118 ]; %rosso scuro, blu scuro, rosso medio, blu medio   
    if export_l == 1
        S.legendl = 0;
    else
        S.legendl = 1;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = 1;
    S.condition_n = 1:4; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
    S.ROIs_labels = {1, 2, 3, 4, 5}; 
    S.lineplot = lineplot;

    S.signtp = [];
    S.signtp_col = [];

    
    S.groups = {1}; % this is one bc you subjects are not divided in groups, set this to 1 and don't ever touch it 
    S.gsubj = {1:83}; %same as above, this containes the indices os subjects and they range from 0 to 83 since your subjects are not divided in groups
    waveplot_groups_local_v2(S) %actual function
    
    if export_l == 1
        exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS2\averagedVSsingle_condition_0' num2str(cc) '.pdf'],'Resolution',300)
    end
    hold on;
end
%% FIGURE S3 VAR EXPLAINED
var_avsub = load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\figureS3\pca_averagedsubjects_COND1_varexp.mat');
var_singlesub = load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS3\pca_singlesubjects_COND1_varexp.mat');
var_concsub = load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS3\pca_concatenatedsubjects_COND1_varexp.mat'); % (?)


var_singlesub = mean(var_singlesub.vare_singlesub(:,:),2);
export_l = 0; % 1 to export the figure
colorline = [0.9 0.4 0.3;0.5 0.4 0.75; 0.3 0.6 0.9];
figure;
plot(var_avsub.vare_avsub(1:20), 'Color',colorline(1,:),'Linewidth', 1.5, 'DisplayName', 'Averaged subjects');
hold on;
plot(var_singlesub(1:20), 'Color',colorline(2,:),'Linewidth', 1.5, 'DisplayName', 'Single subjects');
hold on;
plot(var_concsub.vare(1:20), 'Color',colorline(3,:), 'Linewidth', 1.5 ,'DisplayName', 'Concatenated subjects');
hold on;
grid minor
set(gcf,'color','w')
box on

if export_l == 0
    legend('show');
end

if export_l ==1
    exportgraphics(gcf,'C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS3\Randomization_varexplained.pdf','Resolution',300)
end
%% FIGURE S3: TIME SERIE FROM PCA ON AVERAGED SUBJECTS, BOTH COMPONENTS
% loading data
clear all
addpath('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS3\pca_averagedsubjects_COND1_timeserie.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time

col_l = 0; %1 for significant time-windows with different colors; 1 = only grey color
ylimm = [-1900,800]; %amplitude limits; leave empty [] for automatic adjustment

% conds = 5;  % IN YOUR CASE THESE ARE THE CONDITIONS, NOT 

export_l = 1; %1 = export images; 0 = not

data = permute(J,[2 1 3]);
close all

lineplotdum = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
lineplot = repmat(lineplotdum,1,1);
S = [];
S.sbrim = 0.15;
S.ii = 1;

S.conds = {'1'};
S.data = data(:,:,:); %JUST ONE CONDITION SO NO ITERATIONS
S.STE = 2; %1 = dot lines for standard error; 2 = shadows
S.transp = 0.3; %transparency for standard errors shadow
S.time_real = time_sel(1:1025);
S.colorline = [0.7 0.2 0.2;0.15 0.25 0.6];
if export_l == 1
    S.legendl = 0;
else
    S.legendl = 1;
end
S.x_lim = [-0.1 3.4]; % Set x limits
S.y_lim = ylimm; %Set y limits
S.ROI_n = 1:2;
S.condition_n = 1; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
S.ROIs_labels = {'PC 1' 'PC 1'}; 
S.lineplot = lineplot;

S.signtp = [];
S.signtp_col = [];


S.groups = {1}; % this is one bc you subjects are not divided in groups, set this to 1 and don't ever touch it 
S.gsubj = {1:83}; %same as above, this containes the indices os subjects and they range from 0 to 83 since your subjects are not divided in groups
waveplot_groups_local_v2(S) %actual function
%title(lab(ROIs_to_AAL{ii,1}(pp),:))

if export_l == 1
    exportgraphics(gcf,'C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS3\averagedsubjects_condition_1.pdf','Resolution',300)
end

hold on;
%% FIGURE S3: TIME SERIE FROM PCA ON SINGLE SUBJECTS, BOTH COMPONENTS
% loading data
clear all
addpath('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS3\pca_singlesubjects_timeserie.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time

col_l = 0; %1 for significant time-windows with different colors; 1 = only grey color
ylimm = [-1900,800]; %amplitude limits; leave empty [] for automatic adjustment

% conds = 5;  % IN YOUR CASE THESE ARE THE CONDITIONS, NOT 

export_l = 1; %1 = export images; 0 = not
data = permute(J,[2 1 3]);
close all

lineplotdum = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
lineplot = repmat(lineplotdum,1,1);
S = [];
S.sbrim = 0.15;
S.ii = 1;

S.conds = {'1'};
S.data = data(:,:,:); %JUST ONE CONDITION SO NO ITERATIONS
S.STE = 2; %1 = dot lines for standard error; 2 = shadows
S.transp = 0.3; %transparency for standard errors shadow
S.time_real = time_sel(1:1025);
S.colorline = [0.7 0.2 0.2;0.15 0.25 0.6];          
if export_l == 1
    S.legendl = 0;
else
    S.legendl = 1;
end
S.x_lim = [-0.1 3.4]; % Set x limits
S.y_lim = ylimm; %Set y limits
S.ROI_n = 1:2;
S.condition_n = 1; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
S.ROIs_labels = {'PC 1' 'PC 1'}; 
S.lineplot = lineplot;

S.signtp = [];
S.signtp_col = [];


S.groups = {1}; % this is one bc you subjects are not divided in groups, set this to 1 and don't ever touch it 
S.gsubj = {1:83}; %same as above, this containes the indices os subjects and they range from 0 to 83 since your subjects are not divided in groups
waveplot_groups_local_v2(S) %actual function
if export_l == 1
    exportgraphics(gcf,'C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS3\singlesubjects_condition_1.pdf','Resolution',300)
end

hold on;
%% FIGURE S3 : TIME SERIE OF CONCATENATED SUBJECTS
% loading data
clear all
addpath('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS3\pca_concatenatedsubjects_COND1_timeserie.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time

col_l = 0; %1 for significant time-windows with different colors; 1 = only grey color
ylimm = [-1900,800]; %amplitude limits; leave empty [] for automatic adjustment

% conds = 5;  % IN YOUR CASE THESE ARE THE CONDITIONS, NOT 

export_l = 1; %1 = export images; 0 = not

data = permute(J,[2 1 3]);
close all

lineplotdum = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
lineplot = repmat(lineplotdum,1,1);
S = [];
S.sbrim = 0.15;
S.ii = 1;

S.conds = {'1'};
S.data = data(:,:,:); %JUST ONE CONDITION SO NO ITERATIONS
S.STE = 2; %1 = dot lines for standard error; 2 = shadows
S.transp = 0.3; %transparency for standard errors shadow
S.time_real = time_sel(1:1025);
S.colorline = [0.7 0.2 0.2;0.15 0.25 0.6];        
if export_l == 1
    S.legendl = 0;
else
    S.legendl = 1;
end
S.x_lim = [-0.1 3.4]; % Set x limits
S.y_lim = ylimm; %Set y limits
S.ROI_n = 1:2;
S.condition_n = 1; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
S.ROIs_labels = {'PC 1' 'PC 1'}; 
S.lineplot = lineplot;

S.signtp = [];
S.signtp_col = [];


S.groups = {1}; % this is one bc you subjects are not divided in groups, set this to 1 and don't ever touch it 
S.gsubj = {1:83}; %same as above, this containes the indices os subjects and they range from 0 to 83 since your subjects are not divided in groups
waveplot_groups_local_v2(S) %actual function
%title(lab(ROIs_to_AAL{ii,1}(pp),:))

if export_l == 1
    exportgraphics(gcf,'C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS3\concatenatedsubjects_condition_1.pdf','Resolution',300)
end

hold on;
%% FIGURE S4:
% loading data
addpath('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS4\result01_bonferroni.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS4\result02_fdr.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS4\result03_clusterbased_permutationtest_size.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS4\result04_mcs.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\FigureS4\timeserie_wcoeff_from_maineffect.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time

col_l = 1 ; %1 for significant time-windows with different colors; 0 = only grey color
ylimm = []; %amplitude limits; leave empty [] for automatic adjustment

PCs = 2;  %ROIs (they are 6) IN YOUR CASE THESE ARE THE COMPONENTS, NOT THE REGION OF INTEREST

export_l = 1; %1 = export images; 0 = not
data = permute(J,[2 1 4 3]);
final_stat = [bonferroni_results, fdr_results, Permtest_results, SIGN];
%close all
for cc = 1:PCs %over principal components
    lineplotdum = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
    lineplot = lineplotdum; 
    S = [];
    S.sbrim = 0.15;
    S.ii = cc;

    S.conds = {'Old ','NewT1'};
    S.data = data(cc,:,:,:);
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.time_real = time_sel(1:1025);
    S.colorline = [1 0 0; 0.3686 0.6314 0.7412; 0.8, 0.2, 0.2; 0.6, 0.0, 0.6; 0.0, 0.5, 1.0; 0.0, 0.3, 0.7];        
    if export_l == 1
        S.legendl = 0;
    else
        S.legendl = 1;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = 1;
    S.condition_n = 1:2; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
    S.ROIs_labels = {1, 2}; 
    S.lineplot = lineplot;
    bum = [];
    %%%SIGNIFICATIVITA%%%
    signtp_col = [];
    for ss = 1:4 %over different methods
        sbam = size(final_stat{cc,ss},2);
        clear signtp_temp;
        for ll = 1:sbam %over the number of the significant time-windows for contrast ss
            bum = cat(2,bum,{final_stat{cc,ss}(ll).time_interval}); %concatenating all significant time-windows, for one contrast at a time
            
            signtp_temp(ll) = (ss+2); %getting the color code (number of the significant time-windows for each contrast)
        end
        signtp_col = cat(2, signtp_col, signtp_temp);
    end
    S.signtp = bum;
    if col_l == 1
        S.signtp_col = signtp_col;
        S.colorsign = S.colorline;
    else
        S.signtp_col = [];
        S.colorsign = {'*';'+'};
    end
    
    S.groups = {1}; % this is one bc you subjects are not divided in groups, set this to 1 and don't ever touch it 
    S.gsubj = {1:83}; %same as above, this containes the indices os subjects and they range from 0 to 83 since your subjects are not divided in groups
    waveplot_groups_local_v2(S) %actual function
    
    if export_l == 1
        exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\figureS4\PCAmaineffect_allstats_component_0' num2str(cc) '.pdf'],'Resolution',300)
    end
    hold on;
end

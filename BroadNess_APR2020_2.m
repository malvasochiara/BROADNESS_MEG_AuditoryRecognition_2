%%
%% LBPD functions
% ALWAYS RUN THIS SECTION
%1) Add LBPD functions
% starting up some of the functions for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl); % add to path LBPD and other OSL functions
addpath('/projects/MINDLAB2023_MEG-AuditMemDement/scripts/chiaramalvaso/Broadness')
%% CLUSTER CONFIGURATION
% ALWAYS RUN THIS SECTION IF YOU WANT TO USE CLUSTER
clusterconfig('scheduler', 'cluster');
% clusterconfig('scheduler', 'none'); %use this line to run local
clusterconfig('long_running',1);
clusterconfig('wait',1);
clusterconfig('slot',1); %choose the appropriate amount of memory slot

%% LOADING AVERAGE OVER SUBJECTS
% ALWAYS RUN THIS SECTION BEFORE ANY OTHER SECTION EXCEPT FOR THE SECTIONS: PCA ON SINGLE SUBJECTS AND PCA ON CONCATENATED SUBJECTS

% This section perform the polarity correction on the main effect

%loading data in 3559-voxel space (8mm)
%getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/sources_main_effects.mat');

%actual computation
%adjusting polarity
timex = 45:52;
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if squeeze(mean(t_val_s(jj,timex,1),2)) > 0 %if the data in voxel jj is positive during N100 time
        vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
    else
        vect(jj,1) = 1;
    end
end
dum = zeros(size(t_val_s,1),size(t_val_s,2),size(t_val_s,3));
for cc = 1:size(dum,3) %over conditions
    for jj = 1:size(t_val_s,1) %over brain voxels
        dum(jj,:,cc) = t_val_s(jj,:,cc) .* vect(jj,1); %reversing (or not)..
        disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
    end
end



%% PCA ON MAIN EFFECT (averaged subjects and conditions)
% 1) Time series using PCA on main effect. Saving a matrix with dimensions (time x component x condition x subjects)
% 2) Statistic : table with charactheristics of the significant clusters: 
%       - size
%       - temporal extent
%       - max T value                                                                        
%       - min P value                                                                                            
% 3) Thresholded activation patterns for eachcomponent


%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');

%creating the output directory - choose the approriate name
dirname = 'figure2_TEST';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness';
if ~exist([outputdirpath '/' dirname], 'dir')
    mkdir([outputdirpath '/' dirname])
end
outputdir = [outputdirpath '/' dirname];
%%%%%%%%%%%%%%%%%%%%%%% PCA on main effect to get wcoeff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the necessary input parameters
S = [];
S.H = dum(:,:,1:5); %MAIN EFFECT WITH POLARITY CORRECTION
S.permnum = 1;
S.fig_l = 0; % 1 to plot the results, 0 for no plotting
S.sign_eig = '0';
%S.namenii = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/Test'; %path were you store the .nii images
S.time = time;
S.rand_l = 1;
S.onefig = 0;

[ OUT ] = PCA_LBPD( S );

% extract the values of wcoeffs so you can use them to compute the time series for single subjects
wcoeff = OUT.W;

%%%%%%%%%%%%%%%%%%%%%%%%% Time series for each subject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCs = 1:2; % array containing the principal components you want to consider

list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/SUBJ*.mat');
conds = 5; %number of experimental conditions
J = zeros(size(wcoeff,2),PCs(end),conds,length(list));

for ii = 1:length(list) %over subjects
    load([list(ii).folder '/' list(ii).name])
    t_val_s = OUT.sources_ERFs;
    dum = zeros(size(t_val_s,1),size(t_val_s,2),size(t_val_s,3));
    %adjusting polarity 
    for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(t_val_s,1) %over brain voxels
            dum(jj,:,cc) = t_val_s(jj,:,cc) .* vect(jj,1); %reversing (or not)..
            %            disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
        end
    end
    for iii = 1:conds %over experimental conditions
        J(:,:,iii,ii) = dum(:,1:size(wcoeff,2),iii)' * wcoeff(:,PCs); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
    end
    disp(ii)
end
% saving:
save([outputdir '/timeserie_wcoeff_from_maineffect'], 'J');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Statistic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjust the size of the data and select time points corresponding to the time interval that's relevant for the statistic
starting_time = 0.350; %starting time IN SECONDS
ending_time = 2.50;    %final time IN SECONDS
reduced_index = find(time >= starting_time & time <= ending_time);


dataa = permute(J, [2,1,3,4]);
dataa_lesstime = dataa(:,reduced_index,:,:);
lesstime = time(reduced_index);

comparisons = 4;
PCs = 2;

%check the values to see if they make sense
if comparisons > (size(dataa_lesstime,3)-1)
    error(['ERROR: you want to see ' num2str(comparisons) ' components but you only have ' num2str(size(dataa_lesstime,3)) ' conditions'])
end

if PCs > size(dataa_lesstime,1)
    error(['ERROR: you want to see ' num2str(PCs) ' components but you only computed ' num2str(size(dataa_lesstime,1)) ' components'])
end


% J = (time-points, principal components,conditions,subjects)
P = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)

%%%%%%%%%%% T-test %%%%%%%%%%% 
for ii = 1:size(dataa_lesstime,1) %over principal components
    for jj = 1:size(dataa_lesstime,2) %ove time-points
        for cc = 1:(size(dataa_lesstime,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dataa_lesstime(ii,jj,1,:)),squeeze(dataa_lesstime(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end

end

%%%%%%%%%%% 2)FDR %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear FDR
fdr_siglevel = zeros(PCs,  comparisons);
for comp = 1:PCs %iterate over principal components
    for cond = 1:comparisons % iterate over conditions
        [fdr_siglevel(comp,cond),~,~] = fdr(P(comp,:,cond));
    end
end

fdr_results = cell(PCs,comparisons);

for pp = 1:PCs
    for cc = 1:comparisons
         fdr_matrix = cell(4,length(lesstime)); %initializing the matrix that will be saved on excel
         fdr_matrix(1,:) = num2cell(lesstime);
         clear cluster;
         clear mask;
         mask =double(P(pp,:,cc) < fdr_siglevel(pp,cc));
         cluster = bwconncomp(mask);
         fdr_matrix(2,:) = num2cell(mask);
         for tt = 1:length(lesstime)
             if mask(tt) == 1
                 fdr_matrix{3,tt} = P(pp,tt,cc);
                 fdr_matrix{4,tt} = T(pp,tt,cc);
             else
                 fdr_matrix{3,tt} = 'N.S.';
                 fdr_matrix{4,tt} = 'N.S.';
             end
         end

             filename = [outputdir '/FDR_comp0' num2str(pp) '_comparison_1VS' num2str(cc+1) '.xlsx'];
             writetable(cell2table(fdr_matrix),filename);
         S = {};
         for n = 1:cluster.NumObjects
             idx = cluster.PixelIdxList{n};
             S(n).size = length(idx);
             S(n).time_interval = [lesstime(idx(1)) lesstime(idx(end))]; 
             S(n).min_pval = min(P(pp,idx,cc));
             [absmax, indexmax] = max(abs(T(pp,idx,cc)));
             S(n).Tmax = T(pp,idx(1)-1+indexmax,cc); %storing the Tvalue with the sign
         end
          fdr_results{pp,cc} = S;
    end
end


save([outputdir '/result02_fdr'], 'fdr_results');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thresholded brain images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting weights in the brain (nifti images)
dum_averaged = mean(dum(:,:,1:5),3);
C = cov(dum_averaged');

maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
for ii = 1:PCs(end) %over significant PCs
    wcoeff(:,ii) = OUT.W(:,ii)'*C;
    threshold = mean(abs(wcoeff(:,ii))) + std(abs(wcoeff(:,ii)));
    SO = zeros(size(wcoeff,1),1);
    for jj = 1:size(wcoeff,1)
        if abs(wcoeff(jj,ii)) > threshold
            SO(jj) = wcoeff(jj,ii);
        else
            SO(jj) = 0;
        end
    end
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    for jj = 1:size(SO,1) %over brain sources
        dumm = find(maskk.img == jj); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = SO(jj,:); %storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - component ' num2str(ii)])
    save_nii(nii,[outputdir '/thresholded_actpattern_PCAmaineffect_PC_' num2str(ii) '.nii']); %printing image
end
%% PCA ON RANDOMIZATION OVER TIME
% 1) Performs PCA on randomized data, saving the explained variance and the time series reconstructed on randomized data
% 2) Statistic : table with charactheristics of the significant clusters: 
%       - size
%       - temporal extent                                                                
%       - max T value                                                                        
%       - min P value                                                                  
% 3) Thresholded activation patterns for each component


addpath('/projects/MINDLAB2023_MEG-AuditMemDement/scripts/chiaramalvaso/Broadness');
%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');
%creating the output directory - choose the approriate name
dirname = 'figure3_newcode';
outputdir = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness';
if ~exist([outputdir '/' dirname], 'dir')
    mkdir([outputdir '/' dirname])
end
outputdir = [outputdir '/' dirname];

%%%%%%%%%%%%%%%%% PCA on main effect to get wcoeff for randomization 1: over time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the necessary input parameters
S = [];
S.maineffect = dum(:,:,1:5);
S.folderpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0';
S.outputdir = outputdir; % the function directly store some outputs
S.vect = vect; % copmpute the polarity correcton BEFORE running this section
S.permnum =1000;
S.sign_eig = '0';
S.rand_l = 1; % OVER TIME

job_id = job2cluster(@PCA_RANDOMIZATION, S);
OUT = jobresult(job_id);

% extract the values of wcoeffs to compute the activation patterns
wcoeff = OUT{1,1}.W_RAND;
% extract and save the time series to compute the statistic
J_rand_averaged = OUT{1,1}.timeserie_rand_averaged;
save([outputdir '/randomiazionoverTIME_tsrand_averaged'], 'J_rand_averaged');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Statistic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% demeaning
J_demeaned= zeros(size(J_rand_averaged,1),size(J_rand_averaged,2),size(J_rand_averaged,3),size(J_rand_averaged,4));
for cc = 1:size(J_rand_averaged,4) % over conditions
    for sub = 1:size(J_rand_averaged,3) % over subjects
        for ii = 1:size(J_rand_averaged,2) % over components
            av_signal = mean(J_rand_averaged(:,ii,sub,cc),1);
            J_demeaned(:,ii,sub,cc) = J_rand_averaged(:,ii,sub,cc) - av_signal;
            
        end
    end
end

%adjust the size of the data and select time points corresponding to the desired time interval
starting_time = 0.350; %starting time IN SECONDS
ending_time = 2.50;    %final time IN SECONDS
reduced_index = find(time >= starting_time & time <= ending_time);


dataa = permute(J_demeaned, [2,1,4,3]);
dataa_lesstime = dataa(:,reduced_index,:,:);
lesstime = time(reduced_index);

comparisons = 4;
PCs = 2;

%check the values to see if they make sense
if comparisons > (size(dataa_lesstime,3)-1)
    error(['ERROR: you want to see ' num2str(comparisons) ' comparisons but you only have ' num2str(size(dataa_lesstime,3)) ' conditions'])
end

if PCs > size(dataa_lesstime,1)
    error(['ERROR: you want to see ' num2str(PCs) ' components but you only computed ' num2str(size(dataa_lesstime,1)) ' components'])
end


% J = (time-points, principal components,conditions,subjects)
P = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)

%%%%%%%%%%% T-test %%%%%%%%%%% 
for ii = 1:size(dataa_lesstime,1) %over principal components
    for jj = 1:size(dataa_lesstime,2) %ove time-points
        for cc = 1:(size(dataa_lesstime,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dataa_lesstime(ii,jj,1,:)),squeeze(dataa_lesstime(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end

end

%%%%%%%%%%% 2)FDR %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear FDR
fdr_siglevel = zeros(PCs,  comparisons);
for comp = 1:PCs %iterate over principal components
    for cond = 1:comparisons % iterate over conditions
        [fdr_siglevel(comp,cond),~,~] = fdr(P(comp,:,cond));
    end
end

fdr_results = cell(PCs,comparisons);

for pp = 1:PCs
    for cc = 1:comparisons
         fdr_matrix = cell(4,length(lesstime)); %initializing the matrix that will be saved on excel
         fdr_matrix(1,:) = num2cell(lesstime);
         clear cluster;
         clear mask;
         mask =double(P(pp,:,cc) < fdr_siglevel(pp,cc));
         cluster = bwconncomp(mask);
         fdr_matrix(2,:) = num2cell(mask);
         for tt = 1:length(lesstime)
             if mask(tt) == 1
                 fdr_matrix{3,tt} = P(pp,tt,cc);
                 fdr_matrix{4,tt} = T(pp,tt,cc);
             else
                 fdr_matrix{3,tt} = 'N.S.';
                 fdr_matrix{4,tt} = 'N.S.';
             end
         end

             filename = [outputdir '/randomiazionoverTIME_FDR_comp0' num2str(pp) '_comparison_1VS' num2str(cc+1) '_demeaned.xlsx'];
             writetable(cell2table(fdr_matrix),filename);
         S = {};
         for n = 1:cluster.NumObjects
             idx = cluster.PixelIdxList{n};
             S(n).size = length(idx);
             S(n).time_interval = [lesstime(idx(1)) lesstime(idx(end))]; 
             S(n).min_pval = min(P(pp,idx,cc));
             [absmax, indexmax] = max(abs(T(pp,idx,cc)));
             S(n).Tmax = T(pp,idx(1)-1+indexmax,cc); %storing the Tvalue with the sign
         end
          fdr_results{pp,cc} = S;
    end
end

save([outputdir '/randomiazionoverTIME_result02_fdr_demeaned'], 'fdr_results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thresholded brain images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dum_averaged = mean(dum(:,:,1:5),3);
C = cov(dum_averaged');
%plotting weights in the brain (nifti images)
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
for ii = 1:PCs(end) %over significant PCs
    wcoeff(:,ii) = OUT{1,1}.W_RAND(:,ii)'*C;
    threshold = mean(abs(wcoeff(:,ii))) + std(abs(wcoeff(:,ii)));
    SO = zeros(size(wcoeff,1),1);
    for jj = 1:size(wcoeff,1)
        if abs(wcoeff(jj,ii)) > threshold
            SO(jj) = wcoeff(jj,ii);
        else
            SO(jj) = 0;
        end
    end
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    for jj = 1:size(SO,1) %over brain sources
        dumm = find(maskk.img == jj); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = SO(jj,:); %storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - component ' num2str(ii)])
    save_nii(nii,[outputdir '/randomiazionoverTIME_thresholded_actpattern_PCAmaineffect_PC_' num2str(ii) '.nii']); %printing image
end

%% PCA ON RANDOMIZATION OVER SPACE LABELS
% 1) Performs PCA on randomized data, saving the explained variance and the time series reconstructed on randomized data
% 2) Statistic : table with charactheristics of the significant clusters: 
%       - size
%       - temporal extent  
%       - max T value                                                                        
%       - min P value                                                                                                     
% 3) Thresholded activation patterns for each component


%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');
%creating the output directory - choose the approriate name
dirname = 'figure3_newcode1000perm';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness';
if ~exist([outputdir '/' dirname], 'dir')
    mkdir([outputdir '/' dirname])
end
outputdir = [outputdirpath '/' dirname];
%%%%%%%%%%%%%%%%% PCA on main effect to get wcoeff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = [];
S.maineffect = dum(:,:,1:5);
S.folderpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0';
S.outputdir = outputdir; 
S.vect = vect;
S.permnum =1000;
S.sign_eig = '0';
S.rand_l = 4; % OVER SPACE LABELS

% sending the job to the cluster
job_id = job2cluster(@PCA_RANDOMIZATION, S);
OUT = jobresult(job_id);
% extract the values of wcoeffs to compute the activation patterns
wcoeff = OUT{1,1}.W_RAND;

J_rand = OUT{1,1}.timeserie_rand_averaged;
save([outputdir '/randomiazionoverSPACELABELS_tsrand'], 'J_rand');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Statistic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%adjust the size of the data and select time points corresponding to the desired time interval
starting_time = 0.350; %starting time IN SECONDS
ending_time = 2.50;    %final time IN SECONDS
reduced_index = find(time >= starting_time & time <= ending_time);


dataa = permute(J_rand, [2,1,4,3]);
dataa_lesstime = dataa(:,reduced_index,:,:);
lesstime = time(reduced_index);

comparisons = 4;
PCs = 2;

%check the values to see if they make sense
if comparisons > (size(dataa_lesstime,3)-1)
    error(['ERROR: you want to see ' num2str(comparisons) ' components but you only have ' num2str(size(dataa_lesstime,3)) ' conditions'])
end

if PCs > size(dataa_lesstime,1)
    error(['ERROR: you want to see ' num2str(PCs) ' components but you only computed ' num2str(size(dataa_lesstime,1)) ' components'])
end


% J = (time-points, principal components,conditions,subjects)
P = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)

%%%%%%%%%%% T-test %%%%%%%%%%% 
for ii = 1:size(dataa_lesstime,1) %over principal components
    for jj = 1:size(dataa_lesstime,2) %ove time-points
        for cc = 1:(size(dataa_lesstime,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dataa_lesstime(ii,jj,1,:)),squeeze(dataa_lesstime(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end

end

%%%%%%%%%%% 2)FDR %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear FDR
fdr_siglevel = zeros(PCs,  comparisons);
for comp = 1:PCs %iterate over principal components
    for cond = 1:comparisons % iterate over conditions
        [fdr_siglevel(comp,cond),~,~] = fdr(P(comp,:,cond));
    end
end

fdr_results = cell(PCs,comparisons);

for pp = 1:PCs
    for cc = 1:comparisons
         fdr_matrix = cell(4,length(lesstime)); %initializing the matrix that will be saved on excel
         fdr_matrix(1,:) = num2cell(lesstime);
         clear cluster;
         clear mask;
         mask =double(P(pp,:,cc) < fdr_siglevel(pp,cc));
         cluster = bwconncomp(mask);
         fdr_matrix(2,:) = num2cell(mask);
         for tt = 1:length(lesstime)
             if mask(tt) == 1
                 fdr_matrix{3,tt} = P(pp,tt,cc);
                 fdr_matrix{4,tt} = T(pp,tt,cc);
             else
                 fdr_matrix{3,tt} = 'N.S.';
                 fdr_matrix{4,tt} = 'N.S.';
             end
         end

             filename = [outputdir '/randomiazionoverSPACELABELS_FDR_comp0' num2str(pp) '_comparison_1VS' num2str(cc+1) '.xlsx'];
             writetable(cell2table(fdr_matrix),filename);
         S = {};
         for n = 1:cluster.NumObjects
             idx = cluster.PixelIdxList{n};
             S(n).size = length(idx);
             S(n).time_interval = [lesstime(idx(1)) lesstime(idx(end))]; 
             S(n).min_pval = min(P(pp,idx,cc));
             [absmax, indexmax] = max(abs(T(pp,idx,cc)));
             S(n).Tmax = T(pp,idx(1)-1+indexmax,cc); %storing the Tvalue with the sign
         end
          fdr_results{pp,cc} = S;
    end
end

save([outputdir '/randomiazionoverSPACELABELS_result02_fdr'], 'fdr_results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thresholded brain images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dum_averaged = mean(dum(:,:,1:5),3);
C = cov(dum_averaged');

maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
for ii = 1:PCs(end) %over significant PCs
    wcoeff(:,ii) = OUT{1,1}.W_RAND(:,ii)'*C;
    threshold = mean(abs(wcoeff(:,ii))) + std(abs(wcoeff(:,ii)));
    SO = zeros(size(wcoeff,1),1);
    for jj = 1:size(wcoeff,1)
        if abs(wcoeff(jj,ii)) > threshold
            SO(jj) = wcoeff(jj,ii);
        else
            SO(jj) = 0;
        end
    end
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    for jj = 1:size(SO,1) %over brain sources
        dumm = find(maskk.img == jj); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = SO(jj,:); %storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - component ' num2str(ii)])
    save_nii(nii,[outputdir '/randomiazionoverSPACELABELS_thresholded_actpattern_PCAmaineffectPC_' num2str(ii) '.nii']); %printing image
end


%% PCA ON MAIN EFFECT-saving data for all the significant components
% Same as for the first section PCA ON MAIN EFFECT, but this time a monte carlo simulation is performed
% to get the number of significant principal component and the data are
% stored for all of the significant components


% 1) Time series using PCA on main effect. Saving a matrix with dimensions (time x component x condition x subjects)
% 2) Statistic : table with charactheristics of the significant clusters: 
%       - size
%       - temporal extent             
%       - max T value                                                                        
%       - min P value                                                                                                                                                                                                       
% 3) Thresholded activation patterns for eachcomponent



%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');
%creating the output directory - choose the approriate name
dirname = 'figureS1';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness';
if ~exist([outputdir '/' dirname], 'dir')
    mkdir([outputdir '/' dirname])
end
outputdir = [outputdirpath '/' dirname];
%%%%%%%%%%%%%%%%%%%%%%%PCA on main effect to get wcoeff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the necessary input parameters
S = [];
S.H = dum(:,:,1:5);
S.permnum = 1000; % 1000 permutations in order to see how many PCs we get
S.fig_l = 0; % 1 to plot the results, 0 for no plotting
S.sign_eig = '0';
S.rand_var = 0;
S.rand_ts = 0;
%S.namenii = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/Test'; %path were you store the .nii images
S.time = time;
S.rand_l = 1;
S.onefig = 0;

[ OUT ] = PCA_LBPD( S );

% extract the values of wcoeffs so you can use them to compute the time series for single subjects
wcoeff = OUT.W;
PCs = length(OUT.sign_comps_idx);
disp(['Performing ' num2str(S.permnum) ' permutations, the number of significant components is ' num2str(PCs)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time series for each subject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PCs = 1:length(OUT.sign_comps_idx); % array containing the significant components 
list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/SUBJ*.mat');
conds = 5; % conditions
J = zeros(size(wcoeff,2),PCs(end),conds,length(list));

for ii = 1:length(list) %over subjects
    load([list(ii).folder '/' list(ii).name])
    t_val_s = OUT.sources_ERFs;
    dum = zeros(size(t_val_s,1),size(t_val_s,2),size(t_val_s,3));
    %adjusting polarity 
    for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(t_val_s,1) %over brain voxels
            dum(jj,:,cc) = t_val_s(jj,:,cc) .* vect(jj,1); %reversing (or not)..
            %            disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
        end
    end
    for iii = 1:5 %over experimental conditions
        J(:,:,iii,ii) = dum(:,1:size(wcoeff,2),iii)' * wcoeff(:,PCs); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
    end
    disp(ii)
end
% saving:
save([outputdir '/timeserie_wcoeff_from_maineffect'], 'J');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Statistic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjust the size of the data and select time points corresponding to the desired time interval
starting_time = 0.350; %starting time IN SECONDS
ending_time = 2.50;    %final time IN SECONDS
reduced_index = find(time >= starting_time & time <= ending_time);


dataa = permute(J, [2,1,3,4]);
dataa_lesstime = dataa(:,reduced_index,:,:);
lesstime = time(reduced_index);

comparisons = 4;


%check the values to see if they make sense
if comparisons > (size(dataa_lesstime,3)-1)
    error(['ERROR: you want to see ' num2str(comparisons) ' components but you only have ' num2str(size(dataa_lesstime,3)) ' conditions'])
end

if PCs(end) > size(dataa_lesstime,1)
    error(['ERROR: you want to see ' num2str(PCs(end)) ' components but you only computed ' num2str(size(dataa_lesstime,1)) ' components'])
end


% J = (time-points, principal components,conditions,subjects)
P = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)

%%%%%%%%%%% T-test %%%%%%%%%%% 
for ii = 1:size(dataa_lesstime,1) %over principal components
    for jj = 1:size(dataa_lesstime,2) %ove time-points
        for cc = 1:(size(dataa_lesstime,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dataa_lesstime(ii,jj,1,:)),squeeze(dataa_lesstime(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end

end

%%%%%%%%%%% 2)FDR %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear FDR
fdr_siglevel = zeros(PCs(end),  comparisons);
for comp = 1:PCs(end) %iterate over principal components
    for cond = 1:comparisons % iterate over conditions
        [fdr_siglevel(comp,cond),~,~] = fdr(P(comp,:,cond));
    end
end

fdr_results = cell(PCs(end),comparisons);

for pp = 1:PCs(end)
    for cc = 1:comparisons
         fdr_matrix = cell(4,length(lesstime)); %initializing the matrix that will be saved on excel
         fdr_matrix(1,:) = num2cell(lesstime);
         clear cluster;
         clear mask;
         mask =double(P(pp,:,cc) < fdr_siglevel(pp,cc));
         cluster = bwconncomp(mask);
         fdr_matrix(2,:) = num2cell(mask);
         for tt = 1:length(lesstime)
             if mask(tt) == 1
                 fdr_matrix{3,tt} = P(pp,tt,cc);
                 fdr_matrix{4,tt} = T(pp,tt,cc);
             else
                 fdr_matrix{3,tt} = 'N.S.';
                 fdr_matrix{4,tt} = 'N.S.';
             end
         end

             filename = [outputdir '/FDR_comp0' num2str(pp) '_comparison_1VS' num2str(cc+1) '.xlsx'];
             writetable(cell2table(fdr_matrix),filename);
         S = {};
         for n = 1:cluster.NumObjects
             idx = cluster.PixelIdxList{n};
             S(n).size = length(idx);
             S(n).time_interval = [lesstime(idx(1)) lesstime(idx(end))]; 
             S(n).min_pval = min(P(pp,idx,cc));
             [absmax, indexmax] = max(abs(T(pp,idx,cc)));
             S(n).Tmax = T(pp,idx(1)-1+indexmax,cc); %storing the Tvalue with the sign
         end
          fdr_results{pp,cc} = S;
    end
end


save([outputdir '/result02_fdr'], 'fdr_results');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thresholded brain images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dum_averaged = mean(dum(:,:,1:5),3);
C = cov(dum_averaged');

%plotting weights in the brain (nifti images)
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
for ii = 1:PCs(end) %over significant PCs
    wcoeff(:,ii) = OUT.W(:,ii)'*C;
    threshold = mean(abs(wcoeff(:,ii))) + std(abs(wcoeff(:,ii)));
    SO = zeros(size(wcoeff,1),1);
    for jj = 1:size(wcoeff,1)
        if abs(wcoeff(jj,ii)) > threshold
            SO(jj) = wcoeff(jj,ii);
        else
            SO(jj) = 0;
        end
    end
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    for jj = 1:size(SO,1) %over brain sources
        dumm = find(maskk.img == jj); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = SO(jj,:); %storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - component ' num2str(ii)])
    save_nii(nii,[outputdir '/thresholded_actpattern_PCAmaineffect_PC_' num2str(ii) '.nii']); %printing image
end

%% PCA ON AVERAGED CONDITIONS (and averaged subjects)
% Same as for the first section PCA ON MAIN EFFECT, this section is a reference to compare results with
% the following section, where pca is performed on single conditions

% 1) Time series using PCA on main effect. Saving a matrix with dimensions (time x component x condition x subjects)                                                                        
% 2) Thresholded activation patterns for eachcomponent

%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');
outputdir= '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness';
%creating the output directory - choose the approriate name
dirname = 'figureS2';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness';
if ~exist([outputdir '/' dirname], 'dir')
    mkdir([outputdir '/' dirname])
end
outputdir = [outputdirpath '/' dirname];

% same as before on main effect 
%%%%%%%%%%%%%%%%% PCA on main effect to get wcoeff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the necessary input parameters
S = [];
S.H = dum(:,:,1:5);
S.permnum = 1;
S.fig_l = 0; % 1 to plot the results, 0 for no plotting
S.sign_eig = '0';
S.namenii = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/Test'; %path were you store the .nii images
S.time = time;
S.rand_l = 0; 
S.rand_ts = 0;
S.onefig = 0;

[ OUT ] = PCA_LBPD(S);

% extract the values of wcoeffs so you can use them to compute the time series for single subjects
wcoeff = OUT.W;% extract and saving variance

varexp = OUT.variance_PCS;
save([outputdir '/averagedconditions_varexp_from_maineffect'], 'varexp');

%%%%%%%%%%%%%%%%%%%%%%%%% Time series for each subject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCs = 1:2; % number of components you want to see

list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/SUBJ*.mat');

J = zeros(size(wcoeff,2),PCs(end),5,10);

for ii = 1:length(list) %over subjects
    load([list(ii).folder '/' list(ii).name])
    t_val_s = OUT.sources_ERFs;
    dum = zeros(size(t_val_s,1),size(t_val_s,2),size(t_val_s,3));
    %adjusting polarity 
    for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(t_val_s,1) %over brain voxels
            dum(jj,:,cc) = t_val_s(jj,:,cc) .* vect(jj,1); %reversing (or not)..
            %            disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
        end
    end
    for iii = 1:5 %over experimental conditions
        J(:,:,iii,ii) = dum(:,1:size(wcoeff,2),iii)' * wcoeff(:,PCs); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
    end
    disp(ii)
end
% saving:
save([outputdir '/averagedconditions_timeserie_wcoeff_from_maineffect'], 'J');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thresholded brain images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dum_averaged = mean(dum(:,:,1:5),3);
C = cov(dum_averaged');

maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
for ii = 1:PCs(end) %over significant PCs
    wcoeff(:,ii) = OUT.W(:,ii)'*C;
    threshold = mean(abs(wcoeff(:,ii))) + std(abs(wcoeff(:,ii)));
    SO = zeros(size(wcoeff,1),1);
    for jj = 1:size(wcoeff,1)
        if abs(wcoeff(jj,ii)) > threshold
            SO(jj) = wcoeff(jj,ii);
        else
            SO(jj) = 0;
        end
    end
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    for jj = 1:size(SO,1) %over brain sources
        dumm = find(maskk.img == jj); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = SO(jj,:); %storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - component ' num2str(ii)])
    save_nii(nii,[outputdir '/averagedconditions_thresholded_actpattern_PC_' num2str(ii) '.nii']); %printing image
end
%% PCA ON SINGLE CONDITIONS
% Same as the section before, but this time pca is performed for each
% condition independently

% 1) Time series using PCA on main effect and single condition              
% 2) Thresholded activation patterns for eachcomponent


clear wcoeff;
clear varexp;
%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');
%creating the output directory - choose the approriate name
dirname = 'figureS2';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness';
if ~exist([outputdir '/' dirname], 'dir')
    mkdir([outputdir '/' dirname])
end
outputdir = [outputdirpath '/' dirname];

%%%%%%%%%%%%%%%%% PCA on single conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wcoeff = zeros(size(dum,1),size(dum,2)-1,size(dum,3));
cond = size(dum,3);
varexp = zeros(size(dum,2)-1,size(dum,3));
% Defining the necessary input parameters
S = [];
S.permnum = 1;
S.fig_l = 0; % 1 to plot the results, 0 for no plotting
S.sign_eig = '0';
S.namenii = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/Test'; %path were you store the .nii images
S.time = time;
S.rand_l = 0; % OVER SPATIAL LABELES
S.rand_ts = 0;
S.onefig = 0;
for cc = 1:cond
    S.H = dum(:,:,cc);
    [ OUT ] = PCA_LBPD(S);

    % extract the values of wcoeffs so you can use them to compute the time series for single subjects
    wcoeff(:,:,cc) = OUT.W;% extract and saving variance

    varexp(:,cc) = OUT.variance_PCS;
end
save([outputdir '/singleconditions_varexp_from_maineffect'], 'varexp');

%%%%%%%%%%%%%%%%%%%%%%%%% Time series for each subject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCs = 1:2; %number of components you want to see

list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/SUBJ*.mat');
conds = 5; % number of conditions
J = zeros(size(wcoeff,2),PCs(end),conds,length(list));

for ii = 1:length(list) %over subjects
    load([list(ii).folder '/' list(ii).name])
    t_val_s = OUT.sources_ERFs;
    dum = zeros(size(t_val_s,1),size(t_val_s,2),size(t_val_s,3));
    %adjusting polarity 
    for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(t_val_s,1) %over brain voxels
            dum(jj,:,cc) = t_val_s(jj,:,cc) .* vect(jj,1); %reversing (or not)..
            %            disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
        end
    end
    for iii = 1:5 %over experimental conditions
        J(:,:,iii,ii) = dum(:,1:size(wcoeff,2),iii)' * wcoeff(:,PCs,cc); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
    end
    disp(ii)
end
% saving:
save([outputdir '/singleconditions_timeserie_wcoeff_from_maineffect'], 'J');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thresholded brain images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
for cc = 1:size(dum,3) %over conditions
    clear dum_averaged;
    clear C;
    dum_averaged = mean(dum(:,:,cc),3);
    C = cov(dum_averaged');
    for ii = 1:PCs(end) %over significant PCs
        wcoeff(:,ii,cc) = OUT.W(:,ii,cc)'*C;
        threshold = mean(abs(wcoeff(:,ii,cc))) + std(abs(wcoeff(:,ii,cc)));
        SO = zeros(size(wcoeff,1),1);
        for jj = 1:size(wcoeff,1)
            if abs(wcoeff(jj,ii,cc)) > threshold
                SO(jj) = wcoeff(jj,ii,cc);
            else
                SO(jj) = 0;
            end
        end
        %building nifti image
        SS = size(maskk.img);
        dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
        for jj = 1:size(SO,1) %over brain sources
            dumm = find(maskk.img == jj); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
            [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); %getting subscript in 3D from index
            dumimg(i1,i2,i3,:) = SO(jj,:); %storing values for all time-points in the image matrix
        end
        nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
        nii.img = dumimg; %storing matrix within image structure
        nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
        disp(['saving nifti image - component ' num2str(ii)])
        save_nii(nii,[outputdir '/singleconditions_thresholded_actpattern_condition_' num2str(cc) 'PC_' num2str(ii) '.nii']); %printing image
    end
end
%% PCA ON AVERAGED SUBJECTS and condition = 1
% Same as for the first section PCA ON MAIN EFFECT, this section is a reference to compare results with
% the following two sections, where pca is performed on single subjects and
% concatenated subjects

% 1) Time series using PCA on main effect. Saving a matrix with dimensions (time x component x condition x subjects)       
% 2) Thresholded activation patterns for eachcomponent


%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');
%creating the output directory - choose the approriate name
dirname = 'figureS3';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness';
if ~exist([outputdir '/' dirname], 'dir')
    mkdir([outputdir '/' dirname])
end
outputdir = [outputdirpath '/' dirname];
%%%%%%%%%%%%%%%%%%%%%%%PCA on main effect to get wcoeff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the necessary input parameters

%regularization 
% data_noisy = dum(:,:,1) + 1e-6*rand(size(dum(:,:,1))); % selecting condition 1 and adding a small perturbation
S = [];
S.H =dum(:,:,1); % selecting condition 1
S.permnum = 1;
S.fig_l = 0; % 1 to plot the results, 0 for no plotting
S.sign_eig = '0';
%S.namenii = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/Test'; %path were you store the .nii images
S.rand_var = 0;
S.rand_ts = 0;
S.time = time;
S.rand_l = 1;
S.onefig = 0;

[ OUT ] = PCA_LBPD( S );

% extract the values of wcoeffs so you can use them to compute the time series for single subjects
wcoeff = OUT.W;
vare_avsub = OUT.variance_PCS;
% saving:
save([outputdir '/pca_averagedsubjects_COND1_varexp'], 'vare_avsub');

%%%%%%%%%%%%%%%%%%%%%%%%% Time series for each subject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCs = 1:2; %number of components you want to see

list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/SUBJ*.mat');

J = zeros(size(wcoeff,2),PCs(end),length(list));

for ii = 1:length(list) %over subjects
    load([list(ii).folder '/' list(ii).name])
    t_val_s = OUT.sources_ERFs(:,:,1);
    dum = zeros(size(t_val_s,1),size(t_val_s,2));
    %adjusting polarity 
%     for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(t_val_s,1) %over brain voxels
            dum(jj,:) = t_val_s(jj,:,1) .* vect(jj); %reversing (or not)..
            %            disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
        end
%     end
    J(:,:,ii) = dum(:,1:size(wcoeff,2))' * wcoeff(:,PCs); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
    disp(ii)
end
% saving:
save([outputdir '/pca_averagedsubjects_COND1_timeserie'], 'J');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thresholded brain images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dum_con1 = dum(:,:,1);
C = cov(dum_cond1');

maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
for ii = 1:PCs(end) %over significant PCs
    wcoeff(:,ii) = OUT.W(:,ii)'*C;
    threshold = mean(abs(wcoeff(:,ii))) + std(abs(wcoeff(:,ii)));
    SO = zeros(size(wcoeff,1),1);
    for jj = 1:size(wcoeff,1)
        if abs(wcoeff(jj,ii)) > threshold
            SO(jj) = wcoeff(jj,ii);
        else
            SO(jj) = 0;
        end
    end
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    for jj = 1:size(SO,1) %over brain sources
        dumm = find(maskk.img == jj); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = SO(jj,:); %storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - component ' num2str(ii)])
    save_nii(nii,[outputdir '/pca_averagedsubjects_COND1_thresholded_actpattern_PC_' num2str(ii) '.nii']); %printing image
end

%% PCA ON SINGLE SUBJECTS and condition = 1
% THE THIRD SECTION SHOULD NOT BE RUN BEFORE THIS SECTION
% Same as the section before, but this time pca is computed independently
% for each subject

% 1) Time series using PCA on single subject                                    
% 2) Thresholded activation patterns

%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');
%creating the output directory - choose the approriate name
dirname = 'figureS3';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness';
if ~exist([outputdir '/' dirname], 'dir')
    mkdir([outputdir '/' dirname])
end
outputdir = [outputdirpath '/' dirname];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% adjusting polarity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%loading data in 3559-voxel space (8mm)
%getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/sources_main_effects.mat');
folderpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0';
fileList = dir([folderpath '/SUBJ*.mat']);
N =length(fileList);
bigMatrix = zeros (3559,1026,N); % only consider condition 1
% 
% for i = 1:length(fileList)
for i = 1:N %over subjects
    disp(i);
    currentFileName = fileList(i).name;
    currentFilePath = fullfile(folderpath, currentFileName);
    try
        data = load(currentFilePath);
        bigMatrix(:,:,i) = squeeze(data.OUT.sources_ERFs(:,:,1)); %selcting condition 1
    catch ME
        if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
            disp('errore in %s', currentFilename);
            continue
        end
    end
end
%actual computation
%adjusting polarity: vect is computed on main effect but then the
%correction is applied to each subject independently
timex = 45:52;
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if squeeze(mean(t_val_s(jj,timex,1),2)) > 0 %if the data in voxel jj is positive during N100 time
        vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
    else
        vect(jj,1) = 1;
    end
end

dum = zeros(size(bigMatrix,1),size(bigMatrix,2),N);
for jj = 1:size(t_val_s,1) %over brain voxels
    for s = 1:N % over subjects
        dum(jj,:,s) = bigMatrix(jj,:,s) .* vect(jj,1); %reversing (or not)..
        disp([' - source ' num2str(jj)])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%PCA on main effect to get wcoeff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the necessary input parameters
wcoeff = zeros(size(dum,1),size(dum,2)-1,N);
vare_singlesub = zeros (size(dum,2)-1,N);
S = [];
S.permnum = 1;
S.fig_l = 0; % 1 to plot the results, 0 for no plotting
S.sign_eig = '0';
%S.namenii = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/Test'; %path were you store the .nii images
S.time = time;
S.rand_l = 1;
S.onefig = 0;
for sub = 1:N
    disp(['Performing pca on subject ' num2str(sub)])
    %data_noisy = squeeze(dum(:,:,sub)) + 1e-6*rand(size(squeeze(dum(:,:,sub)))); % selecting condition 1 and adding a small perturbation
    S.H = squeeze(dum(:,:,sub));
    [ OUT ] = PCA_LBPD( S );
    wcoeff(:,:,sub) = OUT.W;
    vare_singlesub(:,sub) = OUT.variance_PCS;
end
save([outputdir '/pca_singlesubjects_COND1_varexp'], 'vare_singlesub');
clear S
clear S_struct
clear OUT
clear bigMatrix
clear vect
clear list
clear timex


%%%%%%%%%%%%%%%%%%%%%%%%% Time series for each subject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCs = 1:2; %number of components you want to see

J = zeros(size(wcoeff,2),PCs(end),N);

for iii = 1:N %over subjects
    J(:,:,iii) = dum(:,1:size(wcoeff,2),iii)' * wcoeff(:,PCs,iii); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
end

% saving:
save([outputdir '/pca_singlesubjects_timeserie'], 'J');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thresholded brain images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dum_averaged = mean(dum(:,:,:),3); % averaging over subjects
C = cov(dum_averaged');

maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure

wcoeff = mean(wcoeff(:,:,:),3); % averaging over subjects
for ii = 1:PCs(end) %over significant PCs
    wcoeff(:,ii) = OUT.W(:,ii)'*C;
    threshold = mean(abs(wcoeff(:,ii))) + std(abs(wcoeff(:,ii)));
    SO = zeros(size(wcoeff,1),1);
    for jj = 1:size(wcoeff,1)
        if abs(wcoeff(jj,ii)) > threshold
            SO(jj) = wcoeff(jj,ii);
        else
            SO(jj) = 0;
        end
    end
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    for jj = 1:size(SO,1) %over brain sources
        dumm = find(maskk.img == jj); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = SO(jj,:); %storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - component ' num2str(ii)])
    save_nii(nii,[outputdir '/pca_singlesubjects_COND1_thresholded_actpattern_PC_' num2str(ii) '.nii']); %printing image
end
%% PCA ON CONCATENATED SUBJECTS and condition = 1
% THE THIRD SECTION SHOULD NOT BE RUN BEFORE THIS SECTION
% Same as the section before, but this time pca is computed on concatenated
% subjects

% 1) Time series using PCA on concatenated subjects                                       
% 2) Thresholded activation patterns



%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');
%creating the output directory - choose the approriate name
dirname = 'figureS3';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness';
if ~exist([outputdir '/' dirname], 'dir')
    mkdir([outputdir '/' dirname])
end
outputdir = [outputdirpath '/' dirname];

%%%%%%%%%%%%%%%%%%%%%%% PCA on concatenated subjects to get wcoeff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIS JOB MUST BE SENT TO THE CLUSTER BECAUSE IT GOES OUT OF MEMORY

% Defining the necessary input parameters
S = [];
S.H = data; 
S.outputdir = outputdir;
S.folderpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0';
S.sign_eig = '0';
jobid = job2cluster(@cluster_concatenated_sub, S);

OUT = jobresult(jobid);
% extract the values of wcoeffs so you can use them to compute the time series for single subjects
wcoeff = OUT{1,1}.W;

clear S
clear OUT

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thresholded brain images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%wcoeff are multiplied by the covariancematrix inside the function @cluster_concatenated_sub
%plotting weights in the brain (nifti images)
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure

wcoeff = mean(wcoeff(:,:,:),3); 
for ii = 1:PCs(end) %over significant PCs
    wcoeff(:,ii) = OUT{1,1}.W(:,ii)'*C;
    threshold = mean(abs(wcoeff(:,ii))) + std(abs(wcoeff(:,ii)));
    SO = zeros(size(wcoeff,1),1);
    for jj = 1:size(wcoeff,1)
        if abs(wcoeff(jj,ii)) > threshold
            SO(jj) = wcoeff(jj,ii);
        else
            SO(jj) = 0;
        end
    end
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    for jj = 1:size(SO,1) %over brain sources
        dumm = find(maskk.img == jj); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = SO(jj,:); %storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - component ' num2str(ii)])
    save_nii(nii,[outputdir '/pca_concatenatedsubjects_COND1_thresholded_actpattern_PC_' num2str(ii) '.nii']); %printing image
end
%% PCA ON MAIN EFFECT - comparison between statistics
% Performs pca on main effect and compute the statistics only for
% comparison condition 1 VS condition 2. Four different types of
% corrections for multiple comparisons are computed: 
%   1) Bonferroni
%   2) FDR
%   3) Cluster based permutation test
%   4) Cluster based MCS
%
% 3) and 4) are computed in the section after this one
%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');
%creating the output directory - choose the approriate name
dirname = 'figureS4';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness';
if ~exist([outputdirpath '/' dirname], 'dir')
    mkdir([outputdirpath '/' dirname])
end
outputdir = [outputdirpath '/' dirname];
%%%%%%%%%%%%%%%%%%%%%%%PCA on main effect to get wcoeff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the necessary input parameters
S = [];
S.H = dum(:,:,1:2); %selcting condition 1 and 2 to compute all the different statistics
S.permnum = 1;
S.fig_l = 0; % 1 to plot the results, 0 for no plotting
S.sign_eig = '0';
%S.namenii = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/Test'; %path were you store the .nii images
S.time = time;
S.rand_l = 1;
S.onefig = 0;

[ OUT ] = PCA_LBPD( S );

% extract the values of wcoeffs so you can use them to compute the time series for single subjects
wcoeff = OUT.W;

%%%%%%%%%%%%%%%%%%%%%%%%% Time series for each subject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCs = 1:2; %number of components you want to see

list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/SUBJ*.mat');

J = zeros(size(wcoeff,2),PCs(end),2, length(list));

for ii = 1:length(list) %over subjects
    load([list(ii).folder '/' list(ii).name])
    t_val_s = OUT.sources_ERFs;
    dum = zeros(size(t_val_s,1),size(t_val_s,2),size(t_val_s,3));
    %adjusting polarity 
    for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(t_val_s,1) %over brain voxels
            dum(jj,:,cc) = t_val_s(jj,:,cc) .* vect(jj,1); %reversing (or not)..
            %            disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
        end
    end
    for iii = 1:2 %over experimental conditions
        J(:,:,iii,ii) = dum(:,1:size(wcoeff,2),iii)' * wcoeff(:,PCs); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
    end
    disp(ii)
end
% saving:
save([outputdir '/timeserie_wcoeff_from_maineffect'], 'J');

%adjust the size of the data and select time points corresponding to the desired time interval
starting_time = 0.350; %starting time IN SECONDS
ending_time = 2.50;    %final time IN SECONDS
reduced_index = find(time >= starting_time & time <= ending_time);

J_lesscond = J(:,:,1:2,:); % if you want to use J computed in section 2
dataa = permute(J_lesscond, [2,1,3,4]);
dataa_lesstime = dataa(:,reduced_index,:,:);
lesstime = time(reduced_index);

comparisons = 1;
PCs = 2;

%check the values to see if they make sense
if comparisons > (size(dataa_lesstime,3)-1)
    error(['ERROR: you want to see ' num2str(comparisons) ' components but you only have ' num2str(size(dataa_lesstime,3)) ' conditions'])
end

if PCs > size(dataa_lesstime,1)
    error(['ERROR: you want to see ' num2str(PCs) ' components but you only computed ' num2str(size(dataa_lesstime,1)) ' components'])
end


% J = (time-points, principal components,conditions,subjects)
P = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)

%%%%%%%%%%% T-test %%%%%%%%%%% 
for ii = 1:size(dataa_lesstime,1) %over principal components
    for jj = 1:size(dataa_lesstime,2) %ove time-points
        for cc = 1:(size(dataa_lesstime,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dataa_lesstime(ii,jj,1,:)),squeeze(dataa_lesstime(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end

end

%%%%%%%%%%% 1)Bonferroni %%%%%%%%%%%

bonferroni_siglevel = 0.05/(size(P,2));

bonferroni_results = cell(PCs,comparisons);

for pp = 1:PCs
    for cc = 1:comparisons
         bonferroni_matrix = cell(4,length(lesstime)); %initializing the matrix that will be saved on excel
         bonferroni_matrix(1,:) = num2cell(lesstime);
         clear cluster;
         clear mask;
         mask = (double(P(pp,:,cc) < bonferroni_siglevel));
         cluster = bwconncomp(mask);
         bonferroni_matrix(2,:) = num2cell(mask);
         for tt = 1:length(lesstime)
             if mask(tt) == 1
                 bonferroni_matrix{3,tt} = P(pp,tt,cc);
                 bonferroni_matrix{4,tt} = T(pp,tt,cc);
             else
                 bonferroni_matrix{3,tt} = 'N.S.';
                 bonferroni_matrix{4,tt} = 'N.S.';
             end
         end

             filename = [outputdir '/Bonferroni_comp0' num2str(pp) '_comparison_1VS' num2str(cc+1) '_finaldata.xlsx'];
             writetable(cell2table(bonferroni_matrix),filename);
         S = {};
         for n = 1:cluster.NumObjects
             idx = cluster.PixelIdxList{n};
             S(n).size = length(idx);
             S(n).time_interval = [lesstime(idx(1)) lesstime(idx(end))]; 
             S(n).min_pval = min(P(pp,idx,cc));
             [absmax, indexmax] = max(abs(T(pp,idx,cc)));
             S(n).Tmax = T(pp,idx(1)-1+indexmax,cc); %storing the Tvalue with the sign
         end
          bonferroni_results{pp,cc} = S;
    end
end

save([outputdir '/result01_bonferroni_finaldata'], 'bonferroni_results');


%%%%%%%%%%% 2)FDR %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear FDR
fdr_siglevel = zeros(PCs,  comparisons);
for comp = 1:PCs %iterate over principal components
    for cond = 1:comparisons % iterate over conditions
        [fdr_siglevel(comp,cond),~,~] = fdr(P(comp,:,cond));
    end
end

fdr_results = cell(PCs,comparisons);

for pp = 1:PCs
    for cc = 1:comparisons
         fdr_matrix = cell(4,length(lesstime)); %initializing the matrix that will be saved on excel
         fdr_matrix(1,:) = num2cell(lesstime);
         clear cluster;
         clear mask;
         mask =double(P(pp,:,cc) < fdr_siglevel(pp,cc));
         cluster = bwconncomp(mask);
         fdr_matrix(2,:) = num2cell(mask);
         for tt = 1:length(lesstime)
             if mask(tt) == 1
                 fdr_matrix{3,tt} = P(pp,tt,cc);
                 fdr_matrix{4,tt} = T(pp,tt,cc);
             else
                 fdr_matrix{3,tt} = 'N.S.';
                 fdr_matrix{4,tt} = 'N.S.';
             end
         end

             filename = [outputdir '/FDR_comp0' num2str(pp) '_comparison_1VS' num2str(cc+1) '_finaldata.xlsx'];
             writetable(cell2table(fdr_matrix),filename);
         S = {};
         for n = 1:cluster.NumObjects
             idx = cluster.PixelIdxList{n};
             S(n).size = length(idx);
             S(n).time_interval = [lesstime(idx(1)) lesstime(idx(end))]; 
             S(n).min_pval = min(P(pp,idx,cc));
             [absmax, indexmax] = max(abs(T(pp,idx,cc)));
             S(n).Tmax = T(pp,idx(1)-1+indexmax,cc); %storing the Tvalue with the sign
         end
          fdr_results{pp,cc} = S;
    end
end


save([outputdir '/result02_fdr_finaldata'], 'fdr_results');


%% CLUSTER- BASED STATISTICS AND CORRECTIONS FOR MULTIPLE COMPARISONS
% Compute the two cluster- based statistics, the previous section MUST be
% run before running this one

addpath('/projects/MINDLAB2023_MEG-AuditMemDement/scripts/chiaramalvaso/Broadness');
% selecting only time points 0.350-2.5 s
%adjust the size of the data and select time points corresponding to the desired time interval
starting_time = 0.350; %starting time IN SECONDS
ending_time = 2.50;    %final time IN SECONDS
reduced_index = find(time >= starting_time & time <= ending_time);


dataa = permute(J_lesscond, [2,1,3,4]); % if you want to use J computed in section 2
dataa_lesstime = dataa(:,reduced_index,:,:);
lesstime = time(reduced_index);

%%%%%%%%%%% T-test %%%%%%%%%%% 
% J = (time-points, principal components,conditions,subjects)
P = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)


for ii = 1:size(dataa_lesstime,1) %over principal components
    for jj = 1:size(dataa_lesstime,2) %ove time-points
        for cc = 1:(size(dataa_lesstime,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dataa_lesstime(ii,jj,1,:)),squeeze(dataa_lesstime(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 3)cluster-based permutation test
comparisons = 1;
PCs = 2;
if PCs > size(dataa_lesstime,1)
    error(['ERROR: you want to see ' num2str(PCs) ' components but you only computed ' num2str(size(dataa_lesstime,1)) ' components'])
end
clear S;
S.nperm  = 1000;
S.alpha = 0.05;
S.threshold = 0.001;
S.time = lesstime;
S.stattype = 'size';
Permtest_results = cell(PCs,comparisons);
elapsedtime_forcomponent = zeros(PCs,1);
for pp = 1:PCs % iterate the process for each principal component
    S.reference = squeeze(dataa_lesstime(pp,:,1,:)); %selecting the current component and condition 1 since we are testing everything against it
    tic;
    for cc = 1:comparisons % iterate the process for each comparisons
        disp(['PC: ' num2str(pp) ' comparison: 1 VS ' num2str(cc+1)])
        
        S.testing = squeeze(dataa_lesstime(pp,:,(cc+1),:)); % selecting the current component and the current condition
        
        Permtest_results{pp,cc} = clusterbased_permutationtest(S); 
    end
    elapsedtime_forcomponent(pp) = toc;
end
save([outputdir '/result03_clusterbased_permutationtest_' (S.stattype) '_finaldata'], 'Permtest_results')
save([outputdir '/result03_clusterbased_permutationtest_' S.stattype '_elapsedtimeforcomponent_finaldata'], 'elapsedtime_forcomponent')



% 4)cluster-based MCS
p_thresh = 0.05; %threshold for binarising p-values vector

SIGN = cell(PCs,comparisons);
elapsedtime_forcomponent = zeros(PCs,1);
for ii = 1:PCs %over components
    tic;
    for cc = 1:comparisons %over conditions 
        Pbin = zeros(1,size(dataa_lesstime,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        tvals = T(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, lesstime, tvals ); %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        %writetable(PDn,[outdir '/AAL_'  AAL_lab{ii} '_OldvsNewT' num2str(cc) '.xlsx'],'Sheet',1); %printing excel file
        SIGN{ii,cc} = sign_clust;
    end
    elapsedtime_forcomponent(ii) = toc;
end

% saving the results in structures just like for the other tests
field_names = {'clustersize','pvalue','time_interval','Tvalue'};

for i = 1:numel(SIGN)
    matrix = SIGN{i};

    S = repmat(struct(),1,size(matrix,1));
    for row = 1:size(matrix,1)
        for col = 1:numel(field_names)
            S(row).(field_names{col}) = matrix{row,col};
        end
    end
    SIGN{i} = S;
end
save([outputdir '/result04_mcs_finaldata'],'SIGN')
save([outputdir '/result04_mcs_elapsedtimeforcomponent_finaldata'], 'elapsedtime_forcomponent')


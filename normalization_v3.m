%% Variables to specify
myFolder = ' '; % Specify the folder where the files are
pattern = ' '; %indicate what the computer should be looking for in the file name, this will define the two cohorts
max_frq = ; %should be either 20 or 60

%% Normalize and plot
% Find files
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return 
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.tsv'); % This can be set to whatever you want, but for power data, keep it on .tsv
theFiles = dir(filePattern); %make a structure that contains all the files for simplicity

% Normalize data
%iterate through every file in the specified folder, sort into control and experimental groups, and do the normalization
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    file = importdata(fullFileName, '\t'); %import tsv as a matrix
    
    if max_frq == 20
        file = file(12:end,5:45); %get rid of useless columns
    elseif max_frq == 60
        file = file(12:end,5:123);
    else
        errorMessage = sprintf('Error: Please define the variable "max_frq" as either 20 or 60');
        uiwait(warndlg(errorMessage));
        break
    end
    
    tf = contains(fullFileName,pattern);

%separate out the different sleep states
    score = file(:,1);
    wake = file(score==1,:);
    nrem = file(score==2,:);
    rem = file(score==3,:);

%normalize data
    if max_frq == 20
        wake_bins = sum(wake(:,2:41));  
        nrem_bins = sum(nrem(:,2:41));
        rem_bins = sum(rem(:,2:41));
    else
        wake_bins = sum(wake(:,2:119));  
        nrem_bins = sum(nrem(:,2:119));
        rem_bins = sum(rem(:,2:119));
    end

    full_wake = sum(wake_bins);
    full_nrem = sum(nrem_bins);
    full_rem = sum(rem_bins);
    
    if tf == 1 %if the expressed pattern is in the file name, save it in one folder, otherwise, save it in another. In this case, the pattern is looking for the control group
        control.normal_wake = (wake_bins./full_wake)*100;
        control.wake(:,:,k) = control.normal_wake; 
        control.normal_nrem = (nrem_bins./full_nrem)*100;
        control.nrem(:,:,k) = control.normal_nrem;
        control.normal_rem = (rem_bins./full_rem)*100;
        control.rem(:,:,k) = control.normal_rem;
    else
        exp.normal_wake = (wake_bins./full_wake)*100;
        exp.wake(:,:,k) = exp.normal_wake;
        exp.normal_nrem = (nrem_bins./full_nrem)*100;
        exp.nrem(:,:,k) = exp.normal_nrem;
        exp.normal_rem = (rem_bins./full_rem)*100;
        exp.rem(:,:,k) = exp.normal_rem; 
    end
end    

% Average processed data
%Control data
control.wake = squeeze(control.wake); %'squeeze' here is turning our 3D matrix into a 2D matrix, which we need for subsequent analysis
control.wake( :, all(~control.wake,1) ) = []; %these will be the files that you want to export to graphpad

control.nrem = squeeze(control.nrem);
control.nrem( :, all(~control.nrem,1) ) = []; %as a side note, what this is doing is taking out all of the empty columns in the matrix

control.rem = squeeze(control.rem);
control.rem( :, all(~control.rem,1) ) = [];

control.avg_wake = mean(control.wake,2); %average all three groups
control.avg_nrem = mean(control.nrem,2);
control.avg_rem = mean(control.rem,2);

%Experimental data
exp.wake = squeeze(exp.wake);
exp.wake( :, all(~exp.wake,1) ) = [];

exp.nrem = squeeze(exp.nrem);
exp.nrem( :, all(~exp.nrem,1) ) = [];

exp.rem = squeeze(exp.rem);
%exp.rem(:,[5 9]) = []; %I'll leave this here as an example for if there are NaN in the data. Just change the numbers in the brackets to which columns contain NaN. Otherwise, comment out
exp.rem( :, all(~exp.rem,1) ) = [];

exp.avg_wake = mean(exp.wake,2);
exp.avg_nrem = mean(exp.nrem,2);
exp.avg_rem = mean(exp.rem,2);

% Area under the curve analysis and statistics
%auc nrem
control.nrem_delta = trapz(control.nrem(1:7,:)); %define frequency interval %Delta- 1:7 Theta- 10:17 Alpha- 18:25 Beta- 26:59 Gamma- 60:118
exp.nrem_delta = trapz(exp.nrem(1:7,:)); %the function 'trapz' is what is actually doing the area under the curve
control.nrem_theta = trapz(control.nrem(10:17,:));
exp.nrem_theta = trapz(exp.nrem(10:17,:)); 
control.nrem_alpha = trapz(control.nrem(18:25,:));
exp.nrem_alpha = trapz(exp.nrem(20:25,:));
control.nrem_beta = trapz(control.nrem(26:59,:));
exp.nrem_beta = trapz(exp.nrem(26:59,:));
control.nrem_gamma = trapz(control.nrem(60:118,:));
exp.nrem_gamma = trapz(exp.nrem(60:118,:));

[stats.h_n_t,stats.p_n_t] = ttest2(control.nrem_theta,exp.nrem_theta); %ttest for significance is default set to .05
[stats.h_n_d,stats.p_n_d] = ttest2(control.nrem_delta,exp.nrem_delta); %the 'h' value is either 1 or 0; 1 if there is significance and 0 if there is not
[stats.h_n_a,stats.p_n_a] = ttest2(control.nrem_alpha,exp.nrem_alpha);
[stats.h_n_b,stats.p_n_b] = ttest2(control.nrem_beta,exp.nrem_beta);
[stats.h_n_g,stats.p_n_g] = ttest2(control.nrem_gamma,exp.nrem_gamma);

%auc wake
control.wake_delta = trapz(control.wake(1:7,:)); 
exp.wake_delta = trapz(exp.wake(1:7,:)); 
control.wake_theta = trapz(control.wake(10:17,:));
exp.wake_theta = trapz(exp.wake(10:17,:));
control.wake_alpha = trapz(control.wake(18:25,:));
exp.wake_alpha = trapz(exp.wake(18:25,:));
control.wake_beta = trapz(control.wake(26:59,:));
exp.wake_beta = trapz(exp.wake(26:59,:));
control.wake_gamma = trapz(control.wake(60:118,:));
exp.wake_gamma = trapz(exp.wake(60:118,:));

[stats.h_w_t,stats.p_w_t] = ttest2(control.wake_theta,exp.wake_theta);
[stats.h_w_d,stats.p_w_d] = ttest2(control.wake_delta,exp.wake_delta);
[stats.h_w_a,stats.p_w_a] = ttest2(control.wake_alpha,exp.wake_alpha);
[stats.h_w_b,stats.p_w_b] = ttest2(control.wake_beta,exp.wake_beta);
[stats.h_w_g,stats.p_w_g] = ttest2(control.wake_gamma,exp.wake_gamma);

%auc rem
control.rem_delta = trapz(control.rem(1:7,:)); 
exp.rem_delta = trapz(exp.rem(1:7,:)); 
control.rem_theta = trapz(control.rem(10:17,:));
exp.rem_theta = trapz(exp.rem(10:17,:)); 
control.rem_alpha= trapz(control.rem(18:25,:));
exp.rem_alpha = trapz(exp.rem(18:25,:));
control.rem_beta = trapz(control.rem(26:59,:));
exp.rem_beta = trapz(exp.rem(26:59,:)); 
control.rem_gamma = trapz(control.rem(60:118,:));
exp.rem_gamma = trapz(exp.rem(60:118,:));

[stats.h_r_t,stats.p_r_t] = ttest2(control.rem_theta,exp.rem_theta);
[stats.h_r_d,stats.p_r_d] = ttest2(control.rem_delta,exp.rem_delta);
[stats.h_r_a,stats.p_r_a] = ttest2(control.rem_alpha,exp.rem_alpha);
[stats.h_r_b,stats.p_r_b] = ttest2(control.rem_beta,exp.rem_beta);
[stats.h_r_g,stats.p_r_g] = ttest2(control.rem_gamma,exp.rem_gamma);

% Plot everything
if max_frq == 20
    bins = (.5:.5:20);
else
    bins = (.5:.5:59);
end

labels={'Delta';'Theta';'Alpha';'Beta';'Gamma'};

figure
if max_frq == 20
    subplot(3,4,1)
    plot(bins,control.avg_wake)
    hold on 
    plot(bins,exp.avg_wake)
    hold off
else
    subplot(3,4,1)
    semilogy(bins,control.avg_wake)
    hold on 
    semilogy(bins,exp.avg_wake)
    hold off
end
ylim([0 10])
legend('control','experimental')
ylabel('Wake % Total Power')
xlabel('Frequency (Hz)')
box off

if max_frq == 20
    subplot(3,4,5)
    plot(bins,control.avg_nrem)
    hold on
    plot(bins,exp.avg_nrem)
    hold off
else
    subplot(3,4,5)
    semilogy(bins,control.avg_nrem)
    hold on
    semilogy(bins,exp.avg_nrem)
    hold off
end
ylim([0 10])
ylabel('NREM % Total Power')
xlabel('Frequency (Hz)')
box off

subplot(3,4,9)
if max_frq == 20
    plot(bins,control.avg_rem)
    hold on
    plot(bins,exp.avg_rem)
    hold off
else
    semilogy(bins,control.avg_rem)
    hold on
    semilogy(bins,exp.avg_rem)
    hold off
end
ylim([0 10])
ylabel('REM % Total Power')
xlabel('Frequency (Hz)')
box off

subplot(3,4,4)
wake_bar = [mean(control.wake_delta),mean(exp.wake_delta);mean(control.wake_theta),mean(exp.wake_theta);mean(control.wake_alpha),mean(exp.wake_alpha);mean(control.wake_beta),mean(exp.wake_beta);mean(control.wake_gamma),mean(exp.wake_gamma)];
bar(wake_bar)
set(gca,'xticklabel',labels)
ylabel('Wake AUC')

subplot(3,4,8)
nrem_bar = [mean(control.nrem_delta),mean(exp.nrem_delta);mean(control.nrem_theta),mean(exp.nrem_theta);mean(control.nrem_alpha),mean(exp.nrem_alpha);mean(control.nrem_beta),mean(exp.nrem_beta);mean(control.nrem_gamma),mean(exp.nrem_gamma)];
bar(nrem_bar)
set(gca,'xticklabel',labels)
ylabel('NREM AUC')

subplot(3,4,12)
rem_bar = [mean(control.rem_delta),mean(exp.rem_delta);mean(control.rem_theta),mean(exp.rem_theta);mean(control.rem_alpha),mean(exp.rem_alpha);mean(control.rem_beta),mean(exp.rem_beta);mean(control.rem_gamma),mean(exp.rem_gamma)];
bar(rem_bar)
set(gca,'xticklabel',labels)
ylabel('REM AUC')

%Make a heatmap
subplot(3,4,2)
imagesc(control.wake)
c = colorbar;
c.Label.String = 'Power (uV^2)';
caxis([0 10])
set(gca,'Ydir','normal')
view([90 -90])
ylabel('Frequency (Hz)')
xlabel('Control wake')
set(gca,'xtick',[])
yt = get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',yt/2)

subplot(3,4,3)
imagesc(exp.wake)
c = colorbar;
c.Label.String = 'Power (uV^2)';
caxis([0 10])
set(gca,'Ydir','normal')
view([90 -90])
ylabel('Frequency (Hz)')
xlabel('Exp wake')
set(gca,'xtick',[])
yt = get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',yt/2)

subplot(3,4,6)
imagesc(control.nrem)
c = colorbar;
c.Label.String = 'Power (uV^2)';
caxis([0 10])
set(gca,'Ydir','normal')
view([90 -90])
ylabel('Frequency (Hz)')
xlabel('Control NREM')
set(gca,'xtick',[])
%colormap(subplot(2,3,2),'winter')
yt = get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',yt/2)

subplot(3,4,7)
imagesc(exp.nrem)
c = colorbar;
c.Label.String = 'Power (uV^2)';
caxis([0 10])
set(gca,'Ydir','normal')
view([90 -90])
ylabel('Frequency (Hz)')
xlabel('Exp NREM')
set(gca,'xtick',[])
%colormap(subplot(2,3,5),'winter')
yt = get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',yt/2)

subplot(3,4,10)
imagesc(control.rem)
c = colorbar;
c.Label.String = 'Power (uV^2)';
caxis([0 10])
set(gca,'Ydir','normal')
view([90 -90])
ylabel('Frequency (Hz)')
xlabel('Control REM')
set(gca,'xtick',[])
%colormap(subplot(2,3,3),'summer')
yt = get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',yt/2)

subplot(3,4,11)
imagesc(exp.rem)
c = colorbar;
c.Label.String = 'Power (uV^2)';
caxis([0 10])
set(gca,'Ydir','normal')
view([90 -90])
ylabel('Frequency (Hz)')
xlabel('Exp REM')
set(gca,'xtick',[])
%colormap(subplot(2,3,6),'summer')
yt = get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',yt/2)
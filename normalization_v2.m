%% Find files
% Specify the folder where the files are
myFolder = 'E:\Projects\04 LC Sleep\LC DREADD\01 analysis\3 hr light\.5 Hz bins light';

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return 
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.tsv'); % This can be set to whatever you want, but for power data, keep it on .tsv
theFiles = dir(filePattern); %make a structure that contains all the files for simplicity

%% Normalize data
%iterate through every file in the specified folder, sort into control and experimental groups, and do the normalization
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    pattern = 'gfp'; %indicate what the computer should be looking for in the file name, this will define the two cohorts
    max_frq = 60; %should be either 20 or 60
    file = importdata(fullFileName, '\t'); %import tsv as a matrix
    
    if max_frq == 20
        file = file(12:end,5:45); %get rid of useless columns
    elseif max_frq == 60
        file = file(12:end,5:90);
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
        wake_bins = sum(wake(:,2:81));  
        nrem_bins = sum(nrem(:,2:81));
        rem_bins = sum(rem(:,2:81));
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
%% Average processed data
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

%% Area under the curve analysis and statistics
%auc nrem
control.nrem_delta = trapz(control.nrem(1:7,:)); %define frequency interval %Delta- 1:7 Theta- 10:19 Alpha- 20:25
exp.nrem_delta = trapz(exp.nrem(1:7,:)); %the function 'trapz' is what is actually doing the area under the curve
control.nrem_theta = trapz(control.nrem(10:19,:));
exp.nrem_theta = trapz(exp.nrem(10:19,:)); 
control.nrem_alpha = trapz(control.nrem(20:25,:));
exp.nrem_alpha = trapz(exp.nrem(20:25,:));

[stats.h_n_t,stats.p_n_t] = ttest2(control.nrem_theta,exp.nrem_theta); %ttest for significance is default set to .05
[stats.h_n_d,stats.p_n_d] = ttest2(control.nrem_delta,exp.nrem_delta); %the 'h' value is either 1 or 0; 1 if there is significance and 0 if there is not
[stats.h_n_a,stats.p_n_a] = ttest2(control.nrem_alpha,exp.nrem_alpha);

%auc wake
control.wake_delta = trapz(control.wake(1:7,:)); 
exp.wake_delta = trapz(exp.wake(1:7,:)); 
control.wake_theta = trapz(control.wake(10:19,:));
exp.wake_theta = trapz(exp.wake(10:19,:));
control.wake_alpha = trapz(control.wake(20:25,:));
exp.wake_alpha = trapz(exp.wake(20:25,:));

[stats.h_w_t,stats.p_w_t] = ttest2(control.wake_theta,exp.wake_theta);
[stats.h_w_d,stats.p_w_d] = ttest2(control.wake_delta,exp.wake_delta);
[stats.h_w_a,stats.p_w_a] = ttest2(control.wake_alpha,exp.wake_alpha);

%auc rem
control.rem_delta = trapz(control.rem(1:7,:)); 
exp.rem_delta = trapz(exp.rem(1:7,:)); 
control.rem_theta = trapz(control.rem(10:19,:));
exp.rem_theta = trapz(exp.rem(10:19,:)); 
control.rem_alpha= trapz(control.rem(20:25,:));
exp.rem_alpha = trapz(exp.rem(20:25,:));

[stats.h_r_t,stats.p_r_t] = ttest2(control.rem_theta,exp.rem_theta);
[stats.h_r_d,stats.p_r_d] = ttest2(control.rem_delta,exp.rem_delta);
[stats.h_r_a,stats.p_r_a] = ttest2(control.rem_alpha,exp.rem_alpha);

%% Plot everything
if max_frq == 20
    bins = (.5:.5:20);
else
    bins = (.5:.5:40);
end

labels={'Delta';'Theta';'Alpha'};

figure
subplot(3,4,1)
plot(bins,control.avg_wake)
hold on 
plot(bins,exp.avg_wake)
hold off
ylim([0 10])
legend('control','experimental')
ylabel('Wake % Total Power')
xlabel('Frequency (Hz)')

subplot(3,4,4)
wake = [mean(control.wake_delta),mean(exp.wake_delta);mean(control.wake_theta),mean(exp.wake_theta);mean(control.wake_alpha),mean(exp.wake_alpha)];
bar(wake)
set(gca,'xticklabel',labels)
ylabel('Wake AUC')

subplot(3,4,5)
plot(bins,control.avg_nrem)
hold on
plot(bins,exp.avg_nrem)
hold off
ylim([0 10])
ylabel('NREM % Total Power')
xlabel('Frequency (Hz)')

subplot(3,4,8)
nrem = [mean(control.nrem_delta),mean(exp.nrem_delta);mean(control.nrem_theta),mean(exp.nrem_theta);mean(control.nrem_alpha),mean(exp.nrem_alpha)];
bar(nrem)
set(gca,'xticklabel',labels)
ylabel('NREM AUC')

subplot(3,4,9)
plot(bins,control.avg_rem)
hold on
plot(bins,exp.avg_rem)
hold off
ylim([0 10])
ylabel('REM % Total Power')
xlabel('Frequency (Hz)')

subplot(3,4,12)
rem = [mean(control.rem_delta),mean(exp.rem_delta);mean(control.rem_theta),mean(exp.rem_theta);mean(control.rem_alpha),mean(exp.rem_alpha)];
bar(rem)
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
% analysis for alternative VPC task with no delay and long delays
% written by Seth Konig 7/2/15. Also for the RM project
%
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Eye Data\';

% datafiles = {'PW150715.2','PW150716.2','PW150717.2'};
% sets = [1 3 5];

% datafiles = {'RR150715.2','RR150716.2','RR150720.2'};
% sets = [1 3 5];
% 
% datafiles = {'RR160627.2'};
% sets = [2];

% datafiles = {'TT150716.2','TT150717.2','TT150720.3'};
% sets = [2 4 6];

% datafiles = {'TO150715.2','TO150716.2','TO150717.2'};
% sets = [2 4 6];



%---Post-lesion---%
% datafiles = {'PW160602.2','PW160603.2','PW160606.2'};
% sets = [2 4 6];

% datafiles = {'RR160627.2','RR160629.2','RR160630.2'};
% sets = [2 4 6];

datafiles = {'TO170822.2','TO170823.2','TO170828.2'};
sets = [1 3 5];

% datafiles = {'PW150715.2','PW150716.2','PW150717.2',...
%     'RR150715.2','RR150716.2','RR150720.2',...
%     'TT150716.2','TT150717.2','TT150720.3',...
%     'TO150715.2','TO150716.2','TO150717.2'};
% sets = [1 3 5 1 3 5 2 4 6 2 4 6];

for file = 1:length(datafiles)
    getVPC_List_EyeData(datafiles{file},sets(file))
end
%%
propleft_propright = cell(1,length(datafiles)); %L/R bias for that session
nov_side = cell(1,length(datafiles)); %the side the repeated or inversely the novel was on
prop_novel = cell(1,length(datafiles));
delay = cell(1,length(datafiles)); %delay between novel and repeat/# of intervening trials
time_nov = cell(1,length(datafiles)); %time course of novelty preference
for file = 1:length(datafiles)
    load([data_dir datafiles{file}(1:8) '_' datafiles{file}(10) '-fixation.mat']);
    
    pl_pr = NaN(2,40);
    del = NaN(1,40);
    n_side = NaN(1,40);
    p_nov = NaN(1,40);
    t_nov = NaN(40,1000);
    for imnum = 1:40
        [row,img_trial] = find(imgnum == imnum);
        if length(row) == 1 %totally novel image during test phase
            
            x = eyedat{img_trial}(1,:);
            y = eyedat{img_trial}(2,:);
            badx = find(x < -16.5 | x > 16.5);
            bady = find(y < -4.5 | y > 4.5);
            
            badind = union(badx,bady);
            x(badind) = [];
            y(badind) = [];
            
            %sometimes cortex lags by a few sample so just shorten to 1000
            if length(x) > 1000
                x = x(1:1000);
                y = y(1:1000);
            end
            
            if row == 1 %then novel image on left
                n_side(imnum) = 1;
                p_nov(imnum) = sum(x < -7.5)./(sum(x < -7.5) + sum(x > 7.5));
                t_nov(imnum,:) = 0;
                t_nov(imnum,x < -7.5) = t_nov(imnum,x < -7.5)+1;
            else %novel image on right
                n_side(imnum) = 2;
                p_nov(imnum) = sum(x > 7.5)./(sum(x < -7.5) + sum(x > 7.5));
                t_nov(imnum,:) = 0;
                t_nov(imnum,x > 7.5) = t_nov(imnum,x > 7.5)+1;
            end
            
            %find delay between fam and test phases
            if row == 1;
                repeat_img = imgnum(2,img_trial);
            else
                repeat_img = imgnum(1,img_trial);
            end
            [row,img_trial] = find(repeat_img == imgnum);
            if length(img_trial) ~= 3
                erorr('Wrong image selected?')
            end
            del(imnum) = img_trial(3)-img_trial(1)-1; %number of intervening trials
            
        else %1st novel side/by side during familiarization phase then repeated during test phase
            img_trial = img_trial(1);
            x = eyedat{img_trial}(1,51:end);
            y = eyedat{img_trial}(2,51:end);
            badx = find(x < -16.5 | x > 16.5);
            bady = find(y < -4.5 | y > 4.5);
            
            badind = union(badx,bady);
            x(badind) = [];
            y(badind) = [];
            
            pl_pr(1,imnum) = sum(x < -7.5);
            pl_pr(2,imnum) = sum(x > 7.5);
            
        end
    end
    
    propleft_propright{file} = pl_pr;
    nov_side{file} = n_side;
    prop_novel{file} =  p_nov;
    delay{file} = del;
    time_nov{file} = t_nov; %time course of novelty preference
end

%calculate the L/R bias
correction = NaN(2,length(datafiles));
for file = 1:length(datafiles);
    mean_pr_pl = nanmean(propleft_propright{file},2);
    correction(1,file) = 2*mean_pr_pl(1)./sum(mean_pr_pl);
    correction(2,file) = 2*mean_pr_pl(2)./sum(mean_pr_pl);
end

correct_time_novel = cell(1,length(datafiles));
correct_prop_novel = cell(1,length(datafiles));
for file = 1:length(datafiles);
    correct_time_novel{file} = NaN(40,1000);
    correct_prop_novel{file} = NaN(1,40);
    for imnum = 1:40
        if ~isnan(nov_side{file}(imnum))
            if nov_side{file}(imnum) ==  1
                correct_prop_novel{file}(imnum)  = prop_novel{file}(imnum)*correction(1,file);
                correct_time_novel{file}(imnum,:) = time_nov{file}(imnum,:)*correction(1,file);
            else
                correct_prop_novel{file}(imnum) = prop_novel{file}(imnum)*correction(2,file);
                correct_time_novel{file}(imnum,:) = time_nov{file}(imnum,:)*correction(2,file);
            end
        end
    end
end

%plot number of L/R novel images
total_left = NaN(1,length(datafiles));
total_right= NaN(1,length(datafiles));
for file = 1:length(datafiles)
    total_left(file) = sum(nov_side{file} == 1);
    total_right(file) = sum(nov_side{file} == 2);
end
figure
bar([total_left' total_right'])
xlabel('Session')
ylabel('Number of Novel Images')
legend('Left','Right')

all_delays = [];
for file = 1:length(datafiles)
    all_delays = [all_delays delay{file}(~isnan(delay{file}))];
end
figure
hist(all_delays,20)
xlabel('# of Intervening stimuli')
ylabel('Count')

mean_pr_pl = NaN(2,length(datafiles));
std_pr_pl =  NaN(2,length(datafiles));
for file = 1:length(datafiles)
    pr_pl = propleft_propright{file};
    pl = pr_pl(1,:);
    pr = pr_pl(2,:);
    pl = 100*pl/sum(nanmean(pr_pl,2));
    pr = 100*pr/sum(nanmean(pr_pl,2));
    mean_pr_pl(1,file) = nanmean(pl);
    mean_pr_pl(2,file) = nanmean(pr);
    std_pr_pl(1,file) = nanstd(pl)/sqrt(20);
    std_pr_pl(2,file) = nanstd(pr)/sqrt(20);
end

figure
hold on
errorb(mean_pr_pl(1,:),std_pr_pl(1,:),'color','blue');
errorb(mean_pr_pl(2,:),std_pr_pl(2,:),'color','red');
plot([0 length(datafiles)+1],[50 50],'k--')
hold off
xlabel('Session')
ylabel('Percent of Time on Each Side')
legend('Left','Right','Chance')
%%
all_time = NaN(length(datafiles),1000);
for file = 1:length(datafiles);
    all_time(file,:) = nanmean(time_nov{file});
end
figure
hold on
dofill(1:5:5000,all_time/10,'blue',1,50)
plot([0 5000],[50 50],'k--')
hold off
xlabel('Time (ms)')
ylabel('Percentage of Time on Novel Image')

%%
delay_index = [0:2:38];
delay_pref = NaN(length(datafiles),20);
for file = 1:length(datafiles)
    for imnum = 1:40;
        if ~isnan(delay{file}(imnum))
            di = find(delay_index == delay{file}(imnum));
            delay_pref(file,di) =  correct_prop_novel{file}(imnum);
        end
    end
end
mean_delay_pref = 100*mean(delay_pref);
std_delay_pref = 100*std(delay_pref)./sqrt(length(datafiles));
figure
hold on
errorb(delay_index,mean_delay_pref,std_delay_pref)
plot([-1 39],[50 50],'k--')
hold off
xlabel('# of intervening stimuli')
ylabel('Novelty Preference')
xlim([-1 39])
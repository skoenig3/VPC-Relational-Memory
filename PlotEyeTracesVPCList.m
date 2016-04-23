data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Eye Data\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\VPClist Figures\';
img_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\VPCListImages\';

datafiles = {'PW150715.2','PW150716.2','PW150717.2',...
    'RR150715.2','RR150716.2','RR150720.2',...
    'TT150716.2','TT150717.2','TT150720.3',...
    'TO150715.2','TO150716.2','TO150717.2'};
sets = [1 3 5 1 3 5 2 4 6 2 4 6];

for file = 6:length(datafiles)
    figure_dir2 = [figure_dir 'Set' num2str(sets(file)) '\'];
    mkdir(figure_dir2)
    img_dir2 =[img_dir 'L' num2str(sets(file)) '\'];
    
    load([data_dir datafiles{file}(1:8) '_' datafiles{file}(10) '-fixation.mat']);
    
    for imnum = 1:20
        [row,img_trial] = find(imgnum == imnum);
        if length(row) == 1 %totally novel image during test phase
            disp('Test  trial,error structure unknown')
            continue
        elseif length(row) ~= 3
            disp('Unknownerror structure unknown')
            continue
        else
            
            nov_img = imgnum(row(1),img_trial(1));
            if row(3) == 1 %novel image is on the opposite row
                rep_img = imgnum(2,img_trial(3));
            else
                rep_img = imgnum(1,img_trial(3));
            end
            img1 = imread([img_dir2  num2str(nov_img) '.bmp']);
            img2 = imread([img_dir2  num2str(rep_img) '.bmp']);
            
            %for familizarition trial
            x = eyedat{img_trial(1)}(1,:);
            y = eyedat{img_trial(1)}(2,:);
            
            badx = find(x < -16.5 | x > 16.5);
            bady = find(y < -4.5 | y > 4.5);
            badind = union(badx,bady);
            x(badind) = [];
            y(badind) = [];
            
            pl = sum(x < -7.5)./(sum(x < -7.5) + sum(x > 7.5));
            pr = sum(x > 7.5)./(sum(x < -7.5) + sum(x > 7.5));
            
            pl_pr = [round(100*pl) round(100*pr)];
            
            x = 24*x+396;
            y = 24*y+108;
            
            img_matrix = 100*ones(216,792,3);
            
            img_matrix(13:204,13:204,:) = img1;
            img_matrix(13:204,589:780,:) = img1;

            subplot(2,1,1)
            imshow(img_matrix/255)
            hold on
            plot(x,y)
            hold off
            axis off
            title(['L: ' num2str(pl_pr(1)) ' R:' num2str(pl_pr(2))])
            
            %for test trial
            x = eyedat{img_trial(3)}(1,:);
            y = eyedat{img_trial(3)}(2,:);
            
            badx = find(x < -16.5 | x > 16.5);
            bady = find(y < -4.5 | y > 4.5);
            badind = union(badx,bady);
            x(badind) = [];
            y(badind) = [];
            
            
            pl = sum(x < -7.5)./(sum(x < -7.5) + sum(x > 7.5));
            pr = sum(x > 7.5)./(sum(x < -7.5) + sum(x > 7.5));
            
            pl_pr = [round(100*pl) round(100*pr)];
            
            x = 24*x+396;
            y = 24*y+108;
            
            img_matrix = 100*ones(216,792,3);
            
            if row(3) == 1 %then novel image on left
                img_matrix(13:204,13:204,:) = img2;
                img_matrix(13:204,589:780,:) = img1;
            else
                img_matrix(13:204,13:204,:) = img1;
                img_matrix(13:204,589:780,:) = img2;
            end

            subplot(2,1,2)
            imshow(img_matrix/255)
            hold on
            plot(x,y)
            hold off
            axis off
            title(['L: ' num2str(pl_pr(1)) ' R:' num2str(pl_pr(2))])
            
            save_and_close_fig(figure_dir2,[datafiles{file}(1:2) '_' num2str(sets(file)) '_VPCList'])     
        end
    end
end
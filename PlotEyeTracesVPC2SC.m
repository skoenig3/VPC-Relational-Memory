%code to plot eye traces on VPC SC images and title with percentage of time
%on image

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Eye Data\';
img_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\VPCSCImages\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Figure\';

% datafiles = {'PW150701.2','PW150706.2','PW150707.2'};
% set =  [1 2 5];
%
% datafiles = {'RR150701.2','RR150706.2','RR150707.2'};
% set = [1 2 5];

% datafiles = {'TO150702.2','TO150706.2','TO150707.2'};
% set =  [3 4 6];

datafiles = {'TT150701.2','TT150702.2','TT150706.2','TT150707.2'};
set = [1 3 2 6];

prop_nov1 = NaN(length(datafiles),20);
prop_nov2 = NaN(length(datafiles),20);
for file = 1:length(datafiles)
    load([data_dir datafiles{file}(1:8) '_' datafiles{file}(10) '-fixation.mat']);
    
    if strcmpi('PW150701.2',datafiles{file})
        eyedat = [eyedat(1:37) {[NaN;NaN]}  eyedat(38:end)];
    elseif strcmpi('TT150701.2',datafiles{file})
        eyedat = [eyedat(1:38) {[NaN;NaN]}  eyedat(39:end)];
    elseif strcmpi('TT150702.2',datafiles{file})
        eyedat = [eyedat(1:38) {[NaN;NaN]}  eyedat(39:end)];
    end
    
    for block = 1:20
        
        figure
        
        for trial = 1:3
            x = eyedat{3*(block-1)+trial}(1,:);
            y = eyedat{3*(block-1)+trial}(2,:);
            
            if isnan(x)
                continue
            end
            
            badx = find(x < -8.75 | x > 8.75);
            bady = find(y < -4.5 | y > 4.5);
            badind = union(badx,bady);
            x(badind) = [];
            y(badind) = [];
            
            pr = sum(x > 0.5);
            pl = sum(x < -0.5);
            
            pl_pr = [round(100*pl/(pr+pl)) round(100*pr/(pr+pl))];
            
            x = 24*x+216;
            y = 24*y+108;
            
            img1 = imread([img_dir 'sc' num2str(set(file)) '\' num2str(imgnum(1,3*(block-1)+trial)) '.bmp']);
            img2 = imread([img_dir 'sc' num2str(set(file)) '\' num2str(imgnum(2,3*(block-1)+trial)) '.bmp']);
            
            img_matrix = 100*ones(216,432,3);
            
            if img_pos(1,3*(block-1)+trial) == -5
                img_matrix(13:204,13:204,:) = img1;
                img_matrix(13:204,229:420,:) = img2;
            else
                img_matrix(13:204,13:204,:) = img2;
                img_matrix(13:204,229:420,:) = img1;
            end
            
            if trial == 2
                if img_pos(1,3*(block-1)+trial) == -5
                    prop_nov1(file,block) =  pl_pr(2);
                else
                    prop_nov1(file,block) =  pl_pr(1);
                end
                
                
            elseif trial == 3
                if img_pos(1,3*(block-1)+trial) == -5
                    prop_nov2(file,block) =  pl_pr(2);
                else
                    prop_nov2(file,block) =  pl_pr(1);
                end
            end
            
            subplot(3,1,trial)
            hold on
            imshow(img_matrix/255)
            plot(x,216-y)
            hold off
            axis off
            title(['%L = ' num2str(pl_pr(1)) '    %R = ' num2str(pl_pr(2))]);
        end
        save_and_close_fig(figure_dir,[datafiles{file}(1:8) 'Block_' num2str(block)])
    end
end

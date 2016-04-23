%set 24 images 53 and 54 are the same
% datafiles = {'PW150609.2','PW150610.2','PW150611.2','PW150612.2',...
%              'PW150615.2','PW150616.2','PW150617.2','PW150618.2',...
%              'PW150619.2','PW150622.2','PW150623.2','PW150624.2',...
%              'PW150625.2','PW150626.2','PW150629.2','PW150630.2'};
% setnums = [41 42 43 44 ...
%            45 46 47 48 ...
%            49 50 51 52 ...
%            53 54 55 56];
% %
% datafiles = {'RR150609.2','RR150610.2','RR150611.2','RR150612.2',...
%              'RR150615.2','RR150616.2','RR150617.2','RR150618.2',...
%              'RR150619.2','RR150622.2','RR150623.2','RR150624.2',...
%              'RR150625.2','RR150626.2','RR150629.2','RR150630.2'};
% setnums = [20 21 22 23 ...
%            24 25 26 27 ...
%            28 29 30 31 ...
%            32 33 34 35];

% datafiles = {'TT150609.2','TT150610.2','TT150611.2','TT150612.2',...
%              'TT150615.2','TT150616.2','TT150617.2','TT150618.2'....
%              'TT150622.2','TT150623.2','TT150624.2','TT150625.2',...
%              'TT150626.2','TT150629.2','TT150630.2'};
% setnums = [20 21 22 23 ...
%            24 25 26 27 ...
%            28 29 30 31 ...
%            32 33 34];

datafiles = {'TO150609.2','TO150610.2','TO150611.2','TO150612.2',...
    'TO150615.2','TO150616.2','TO150617.2','TO150618.2'....
    'TO150622.2','TO150623.2','TO150624.2','TO150625.2',...
    'TO150626.2','TO150629.2','TO150630.2'};
setnums = [20 21 22 23 ...
    24 25 26 27 ...
    28 29 30 31 ...
    32 33 34];

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Eye Data\';
img_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\VPCImages\';

clr = ['rgbmk'];
for file = 1:length(datafiles)
    load([data_dir datafiles{file}(1:8) '_' datafiles{file}(10) '-fixation.mat']);
    
    for block = 1:20
        
        figure
        
        for trial = 1:2
            x = eyedat{2*(block-1)+trial}(1,:);
            y = eyedat{2*(block-1)+trial}(2,:);
            
            badx = find(x < -8.75 | x > 8.75);
            bady = find(y < -4.5 | y > 4.5);
            
            badind = union(badx,bady);
            x(badind) = [];
            y(badind) = [];
            
            x = x(1:1000);
            y = y(1:1000);
            
            x = 24*x+216;
            y = 24*y+108;
            
            img1 = imread([img_dir '2im0' num2str(setnums(file)) '\' num2str(imgnum(1,2*(block-1)+trial)) '.bmp']);
            img2 = imread([img_dir '2im0' num2str(setnums(file)) '\' num2str(imgnum(2,2*(block-1)+trial)) '.bmp']);
            
            img_matrix = 100*ones(216,432,3);
            
            if img_pos(1,trial) == -5
                img_matrix(13:204,13:204,:) = img1;
                img_matrix(13:204,229:420,:) = img2;
            else
                img_matrix(13:204,13:204,:) = img2;
                img_matrix(13:204,229:420,:) = img1;
            end
            
            subplot(2,1,trial)
            imshow(img_matrix/255)
            hold on
            for t = 1:5
                tim = 200*(t-1)+1 : t*200; 
                plot(x(tim),216-y(tim),clr(t))
            end
            hold off
            axis off
            
        end
    end
end
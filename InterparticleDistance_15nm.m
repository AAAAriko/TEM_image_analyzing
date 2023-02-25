%% 15nm-SiO2 particle size distribution analysis
% by Jirameth 'Mett' Tarnsangpradit
% last update: 2/24/2023 
%
% To determine the size distribution of small particles from TEM images

clf
close all
clear

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Loading image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[file,path] = uigetfile({'*.jpg;*.fig;*.png;*.tif;*'}); %load file
 disp(['File: ',path,file])
[pathstr,name,ext] = fileparts(file); 
s1 = name;
s2 = ext;
a0 = imread(strcat(path,s1,s2));
if size(a0,3) >= 3 %checking if the loaded image has 3 (RGB) or 1 channel
    a0 = rgb2gray(a0); %convert to single channel
end
figure(1)
set(gcf, 'Units','normalized','Position',[0.01,0.5,0.4,0.4]) % set where the figure should pop up on the screen
subplot(6,9,[1:3,10:12,19:21]) %first 2 number is figure size, the rest is the array where the plot should lie in figure
imshow(a0)
title('Figure1.1: Original Image');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Image cleaning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row = size(a0,1);
col = size(a0,2);
sigma = floor(0.1/100*row); % determines how much the image is blurred
a1 = 255-imgaussfilt(255-a0,sigma);

% subplot(1,2,2)
subplot(6,9,[4:6,13:15,22:24])
imshow(a1)
title('Figure1.2: Smoothed Image');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Scale bar check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(6,9,[7:9,16:18])
bar_sens = 0.5;
a2 = imbinarize(255-a0, 'adaptive', 'ForegroundPolarity', 'dark', 'Sensitivity',bar_sens);
a2_1 = 1-a2;
imshow(a2_1(row*3/4:row,1:col*2/4))
title('Figure2.1: Binarize')

subplot(6,9,[25:27,34,36])
a2_open_sens = round(row*col*(1e-3));
a2_2 = bwareaopen(a2_1,a2_open_sens); % removing white area with less than # of specified pixels
imshow(a2_2(row*3/4:row,1:col*2/4))
title('Figure2.2: Open')

subplot(6,9,[43:45,52:54])
se_a2 = strel('square', 25);
a2_3 = imerode(a2_2, se_a2);  % erode to remove refine edge
a2_3 = imdilate(a2_3, se_a2); % dilate to restore removed edge
imshow(a2_3(row*3/4:row,1:col*2/4))
title('Figure2.3: Erode & dilate')

[B,L,n] = bwboundaries(a2_3); % scan for object & find size
scale_px = max(B{1}(:,2)) - min(B{1}(:,2))+1; % find scale bar length in px

% dialog box for scale size in nm
prompt = {'Enter Scale bar size (nm):'};
dlgtitle = 'Input';
dims = [1,35]; % dialog box dimension
definput = {''}; %default value in dialog box
scale_nm = inputdlg(prompt,dlgtitle,dims,definput);
scale_nm = str2double(scale_nm{1});
nmperpx = scale_nm/scale_px; %determine length in nm per pixel
disp(['Number of pixel = ', num2str(scale_px)])
disp(['distance per pixel = ',num2str(nmperpx),' nm/pixel'])
%}

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Particles detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4.1. Image processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(6,9,[28:30,37:39,46:48])
a3_1 = imgaussfilt(255-a0,8);
a3_1 = imbinarize(255-a3_1, 'adaptive', 'ForegroundPolarity', 'dark', 'Sensitivity',0.5);
imshow(1-a3_1)
title('Figure3.1: Binarize')

subplot(6,9,[31:33,40:42,49:51])
a3_2 = bwareaopen(1-a3_1,2000);
a3_2 = imfill(a3_2, 'holes');
imshow(a3_2)
title('Figure3.2: Open')

% 4.2. Objects labeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[B,L,n,A] = bwboundaries(a3_2); % find the boundaries of the particles

figure(2)
imshow(a3_2)
set(gcf, 'Units','normalized','Position',[0.42,0.05,0.5,0.85])
hold on

colors=['b' 'g' 'r' 'c' 'm' 'y'];
for k = 1:length(B)
    boundary = B{k};
    cidx = mod(k,length(colors))+1;
    plot(boundary(:,2),boundary(:,1),colors(cidx),'LineWidth',1);
    visboundaries(B{k})
    rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
    col = boundary(rndRow,2);
    row = boundary(rndRow,1);
    h = text(col,row,num2str(L(row,col)));
    set(h,'Color',colors(cidx),'FontSize',16);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Distribution calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5.1. Obtaining particle sizes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The goal here is to remove objects that are not particle or 
% overlappping particles that are not distinguishable via image
% processing in the previous steps

stats = regionprops(a3_2,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity');
% while loop to repeatedly removing non-particle object
loop = 1;
while loop > 0
    if loop == 1
        % dialog box for particle removal
        prompt = {'Enter particle(s) number to be removed:'};
        dlgtitle = 'Input';
        dims = [1,40]; % dialog box dimension
        definput = {''}; %default value in dialog box
        non_particle = inputdlg(prompt,dlgtitle,dims,definput);
        non_particle = str2num(convertCharsToStrings(non_particle{1})); %#ok<ST2NM>
        del = non_particle;
        disp(['Removed particles:',num2str(del)])
        
        % show original image
        close(figure(2))
        figure(3)
        set(gcf, 'Units','normalized','Position',[0.01,0.05,0.95,0.85])
        subplot(1,3,1)
        imshow(a1)
        title('Gaussian blur')
        % show pre-object-removal image
        subplot(1,3,2)
        imshow(a3_2)
        hold on
        for k = 1:length(B)
            boundary = B{k};
            cidx = mod(k,length(colors))+1;
            plot(boundary(:,2),boundary(:,1),colors(cidx),'LineWidth',1);
            visboundaries(B{k})
            rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
            col = boundary(rndRow,2);
            row = boundary(rndRow,1);
            h = text(col,row,num2str(L(row,col)));
            set(h,'Color',colors(cidx),'FontSize',12);
        end
        title('Pre-object-removal')
        % show post-object-removal image
        subplot(1,3,3)
        imshow(a3_2)
        hold on
        for k = 1:length(B)
            if k~= del
                boundary = B{k};
                cidx = mod(k,length(colors))+1;
                plot(boundary(:,2),boundary(:,1),colors(cidx),'LineWidth',1);
                visboundaries(B{k})
                rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
                col = boundary(rndRow,2);
                row = boundary(rndRow,1);
                h = text(col,row,num2str(L(row,col)));
                set(h,'Color',colors(cidx),'FontSize',12);
            end
        end
        title('Post-object-removal')
        loop = loop + 1;
    else
        % 
        answer = questdlg('Any other particle(s) should be removed?','','Remove','Add','No','Remove');
        switch answer
            case 'Remove'
                % prompt dialog box and add to removal list
                prompt = {'Enter particle(s) number to be removed:'};
                dlgtitle = 'Input';
                dims = [1,40]; % dialog box dimension
                definput = {''}; %default value in dialog box
                non_particle = inputdlg(prompt,dlgtitle,dims,definput);
                non_particle = str2num(convertCharsToStrings(non_particle{1})); %#ok<ST2NM>
                if ~isempty(non_particle)
                    del = union(del,non_particle);
                    disp(['Removed particles:',num2str(del)])
                end
                % replot figure3.3
                subplot(1,3,3)
                imshow(a3_2)
                for k = 1:length(B)
                    if k~= del
                        boundary = B{k};
                        cidx = mod(k,length(colors))+1;
                        plot(boundary(:,2),boundary(:,1),colors(cidx),'LineWidth',1);
                        visboundaries(B{k})
                        rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
                        col = boundary(rndRow,2);
                        row = boundary(rndRow,1);
                        h = text(col,row,num2str(L(row,col)));
                        set(h,'Color',colors(cidx),'FontSize',12);
                    end
                end
            case 'Add'
                prompt = {'Enter particle(s) number to be removed:'};
                dlgtitle = 'Input';
                dims = [1,40]; % dialog box dimension
                definput = {''}; %default value in dialog box
                add_particle = inputdlg(prompt,dlgtitle,dims,definput);
                add_particle = str2num(convertCharsToStrings(add_particle{1})); %#ok<ST2NM>
                if ~isempty(add_particle)
                    del = setdiff(del,add_particle);
                    disp(['Removed particles:',num2str(del)])
                end
                % replot figure3.3
                subplot(1,3,3)
                imshow(a3_2)
                for k = 1:length(B)
                    if k~= del
                        boundary = B{k};
                        cidx = mod(k,length(colors))+1;
                        plot(boundary(:,2),boundary(:,1),colors(cidx),'LineWidth',1);
                        visboundaries(B{k})
                        rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
                        col = boundary(rndRow,2);
                        row = boundary(rndRow,1);
                        h = text(col,row,num2str(L(row,col)));
                        set(h,'Color',colors(cidx),'FontSize',12);
                    end
                end
%                 loop = 0;
                
            case 'No'
                loop = 0;
        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Find interparticle distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The interparticle distance can be calculated using the centroid
% information obtained from 'regionprops' above. We then determine the
% nearest neighbor and calculate their distance. Once center-to-center
% distance is calculated, we can subtract the value by each of the radius
% and plot the distribution.

%{
Ctrd = vertcat(stats.Centroid);
D = sqrt((Ctrd(:,1)-Ctrd(:,1).').^2 + (Ctrd(:,2)-Ctrd(:,2).').^2);
D(abs(D)<1E-6) = NaN;
[Dn Nn] = min(D,[],2); 
% Dn is distance to nearest neighbor
% Nn is nearest neighbor
% plot lines between points and the calculated nearest neighbor


% loop over D to find:
% 1) k-Nearest neighbor distances (mink?)
% 2) index of all k neighbors
% Note: check for non-particle and skip over that

for i
%}


figure(4)
imshow(a3_2)
hold on

remain_B = B;
remain_B(del) = [];
for k = 1:length(remain_B)
    boundary = remain_B{k};
    cidx = mod(k,length(colors))+1;
    plot(boundary(:,2),boundary(:,1),colors(cidx),'LineWidth',1);
    visboundaries(remain_B{k})
    rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
    col = boundary(rndRow,2);
    row = boundary(rndRow,1);
    h = text(col,row,num2str(L(row,col)));
    set(h,'Color',colors(cidx),'FontSize',12);
end

%{
remain_Ctrd = Ctrd;
remain_Ctrd(del,:) = [];

for i = 1:size(Ctrd,1)
    if 
    [Idx,D] = knnsearch(Ctrd([1:i-1,i+1:end],:),Ctrd);
    end
end
%}


% references https://people.sc.fsu.edu/~jburkardt/presentations/voronoi_neighbors.pdf
% Determine voronoi cell nearest neighbor (show nxn array)
Ctrd = vertcat(stats.Centroid);
remain_Ctrd = Ctrd;
remain_Ctrd(del,:) = [];
[V,C] = voronoin(remain_Ctrd);

% find edges & find cell that share such edge
% 2 particles are neighbor if they share and edge or 2 vertices
n = length(C);
% vn = sparse(n,n);
vn = zeros(n,n);
interdist = [];
for i = 1:n
    for j = i+1:n
        s = size(intersect(C{i},C{j})); 
        if s(2) > 1
            vn(j,i) = 1;
            vn(i,j) = 1;
            interdist(end+1) = sqrt((remain_Ctrd(i,1)-remain_Ctrd(j,1).').^2 + (remain_Ctrd(i,2)-remain_Ctrd(j,2).').^2);
        end 
    end
end

voronoi(remain_Ctrd(:,1),remain_Ctrd(:,2),'-y')

interdist = nmperpx*interdist.';

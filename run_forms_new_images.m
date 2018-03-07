%% EXAMPLE SCRIPT FOR FITTING A 3D MORPHABLE MODEL TO A SET OF IMAGES
%
% Jason Manley, updated Mar 2018
% Question? jmanley@rockefeller.edu

addpath(genpath('/path/to/repo/forms/'))

path = '/path/to/your/images/';
file_type = '.your_file_type';

files = dir(fullfile(path,['*' file_type]));
N = length(files);

project = struct();


%% load your own mesh variables, or use one of those provided by Cashman & Fitzgibbon

% example: use their dolphins mesh
dolphins = read_project('projects/dolphins.fpj');

project.vertices    = dolphins.vertices;
project.mesh        = dolphins.mesh;
project.faces       = dolphins.faces;
project.cand_ixs    = dolphins.cand_ixs;
project.cand_uvs    = dolphins.cand_uvs;
project.cand_limits = dolphins.cand_limits;
project.cand_dists  = dolphins.cand_dists;
project.cand_derivs = dolphins.cand_derivs;


%% load images

project.images = struct();

for i=1:length(files)
    project.images(i).image = imread(files(i).name);
end


%% find silhouettes
%  here you should use a method for finding the silhouette of your object
%  of interest.

% simple example: convert images to grayscale, threshold images, find 
% silhouette around largest region

threshold = 20; % this will need to be tweaked for your images...

figure;
for i=1:N
    grayscale = rgb2gray(project.images(i).image);
    props = regionprops(grayscale>threshold,'PixelIdxList','Area');
    [~,sortIdx] = sort(props.Area,'descend');
    
    % find largest region
    mask = zeros(size(grayscale));
    mask(props(sortIdx(1)).PixelIdxList) = 1; 

    % find silhouette
    [i,j] = ind2sub(size(mask),props(sortIdx(1)).PixelIdxList(1));
    project.images(i).silhouette = bwtraceboundary(mask,[i j], 'E');
    
    % plot result
    clf; 
    imagesc(project.images(i).image); title(['Image ' num2str(i)])
    plot(project.images(i).silhouette(:,2),images(i).silhouette(:,1),'r','LineWidth',2);
    pause;
end


%% set initial conditions for all images

for i=1:N
    dolphins = set_initial_conditions(project,i);
end

% if curious, plot all dolphins with silhouettes and constraint points
f = figure;
for i=1:N
    clf;
    imshow(project.images(i).image,'Border','tight'); title(['Image ' num2str(i)])
    hold on; plot(project.images(i).silhouette(:,1),project.images(i).silhouette(:,2),'r','linewidth',3)
    scatter(project.images(i).constraints2d(:,1),project.images(i).constraints2d(:,2),'w','filled')
end

% if curious, plot initial conditions
figure; 
for i=1:N
    clf;
    plot_modelfit(project,i); title(['Image ' num2str(i)])
    pause;
end


%% fit the model!

parameters = makeModelParameters; % likely should test various parameter values!

model = forms(project, parameters);


%% visualize model fits

% model fits by image
figure;
for i=1:N
     clf;
    subplot(1,3,1); imagesc(project.images(i).image);
    hold on; scatter(project.images(i).silhouette(:,1),project.images(i).silhouette(:,2),10,'r','filled')
    title('Original Frame','fontsize',14); axis off;
    
    a=subplot(1,3,2); plot_modelfit(project,i,[],true,false); axis off;
    title('Initial Conditions','fontsize',14);
    
    b=subplot(1,3,3); plot_modelfit(project,i,dolphin_model,true,false); axis off;
    title('Model Fit','fontsize',14)
    linkaxes([a,b],'xy');
    Link = linkprop([a, b], ...
       {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
    setappdata(gcf, 'StoreTheLink', Link);
    pause;
end

% basis shapes
figure;
for i=1:parameters.M
    plot_basisshapenorm(dolphins,dolphin_model,i); title(['Basis Shape ' num2str(i)])
    pause;
end




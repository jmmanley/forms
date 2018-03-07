function transformed = plot_modelfit_video(project, i, model, videoPath, showsilhouette, plotmesh)

% PLOT_MODELFIT  Plots the fit of the model for instance 'i'.
%   
%   To view the initial fit before optimization, do not pass 'model'.
%
% Makes a video where the camera revolves around the 3D fit.
%
% Jason Manley, Nov 2017

if nargin < 5
    showsilhouette = true;
end

if nargin < 6
    plotmesh = false;
end

f = figure('visible','off','position',[1 1 1000 300]);
subplot(1,3,1);
imshow(project.images(i).image);
hold on; plot(project.images(i).silhouette(:,1),project.images(i).silhouette(:,2),'r','linewidth',2);
title('Original Frame','fontsize',14); axis off;

a=subplot(1,3,2);
plot_modelfit(project,i,[],showsilhouette,plotmesh);
title('Initial Conditions','fontsize',14); axis vis3d off;

b=subplot(1,3,3);
plot_modelfit(project,i,model,showsilhouette,plotmesh);
 title('Model Fit','fontsize',14); axis vis3d off;

linkaxes([a,b],'xy');
Link = linkprop([a, b], ...
       {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
    setappdata(f, 'StoreTheLink', Link);
    
vw = VideoWriter([videoPath 'fit_' num2str(i) '.mp4'],'MPEG-4');
vw.open;

dtheta = 4;
loops = 5;
numFrames = 360/dtheta * loops;

for frame=1:numFrames
    writeVideo(vw,getframe(f));
    camorbit(a,dtheta,0,'data',[0 1 0]);
end

vw.close;
close(f);


end


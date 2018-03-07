function [project] = set_initial_conditions(project,i,onlyConstraints)
% A tool for inputing the initial conditions for the Cashman & Fitzgibbon,
% 2012 method from "What Shape Are Dolphins? Building 3D Morphable Models
% from 2D Images"
% Jason Manley, updated Oct 2017
% Questions? jmanley@rockefeller.edu

if nargin < 3
    onlyConstraints = false;
end


% figure settings
img_width = 400;
img_height = 330;
img_sep = 50;
fig_margin = 20;
slide_width = 120;
slide_sep = 20;
slide_sep_height = 40;


%% Initialization

% create window
fig_width = 2*img_width+img_sep+fig_margin+slide_sep+slide_width;
fig_height = img_height+fig_margin;
f = figure('visible','off','Resize','off','position', [1000,1000,fig_width,fig_height],...
    'NumberTitle','off','name','Set Initial Conditions','WindowKeyPressFcn',...
    @keypress_callback);

% left image
l_ax = axes('Units','Pixels','XTick',[],'YTick',[],...
    'Position', [fig_margin, fig_margin, img_width, img_height]);
l_img = imshow(project.images(i).image); hold on; scatter(project.images(i).silhouette(:,1),project.images(i).silhouette(:,2),15,'r')
set(l_ax,'ydir','reverse');

% right image
r_ax = axes('Units','Pixels','XTick',[],'YTick',[],...
    'Position', [fig_margin+img_width+img_sep, fig_margin, img_width, img_height]);
plot_looplimit(project.mesh,project.vertices); axis off
hold on; scatter3(project.vertices(:,1),project.vertices(:,2),project.vertices(:,3),'filled')
set(r_ax,'ydir','reverse')
x=0;
y=0;
z=0;

sld_x = uicontrol('Style', 'slider',...
        'Min',-180,'Max',180,'Value',0,...
        'Position', [fig_margin+2*img_width+img_sep+slide_sep img_height slide_width 20],...
        'Callback', @update_x); 
    
sld_y = uicontrol('Style', 'slider',...
        'Min',-180,'Max',180,'Value',0,...
        'Position', [fig_margin+2*img_width+img_sep+slide_sep img_height-slide_sep_height slide_width 20],...
        'Callback', @update_y); 
    
sld_z = uicontrol('Style', 'slider',...
	'Min',-180,'Max',180,'Value',0,...
	'Position', [fig_margin+2*img_width+img_sep+slide_sep img_height-slide_sep_height*2 slide_width 20],...
	'Callback', @update_z);
    
    function update_x(source,event)
        x = source.Value*pi/180;
        rotateCamera;
    end

    function update_y(source,event)
        y = source.Value*pi/180;
        rotateCamera;
    end

    function update_z(source,event)
        z = source.Value*pi/180;
        rotateCamera;
    end
    
    % plot from new camera angle
    function rotateCamera
        R = rotMatrix(x,y,z);
        newVerts = (R(1:3,1:3)*project.vertices')';
        cla(r_ax); plot_looplimit(project.mesh,newVerts); axis off
        hold on; scatter3(newVerts(:,1),newVerts(:,2),newVerts(:,3),'filled')
        scatter3(newVerts(86,1),newVerts(86,2),newVerts(86,3),'r','filled')
        scatter3(newVerts(87,1),newVerts(87,2),newVerts(87,3),'k','filled')
        scatter3(newVerts(99,1),newVerts(99,2),newVerts(99,3),'r','filled')
        scatter3(newVerts(100,1),newVerts(100,2),newVerts(100,3),'k','filled')
        scatter3(newVerts(42,1),newVerts(42,2),newVerts(42,3),'c','filled')
        scatter3(newVerts(90,1),newVerts(90,2),newVerts(90,3),'c','filled')
        set(r_ax,'ydir','reverse')
        drawnow;
    end

rotateCamera;
set(f,'Visible','On')


%% Manually rotate to find azimuth and elevation

if ~onlyConstraints
    disp('Rotate dolphin mesh to match the approximate position of the image on the left.')
    input('Press "Enter" to continue...')
end


%% Find constraint points

keepGoing = true;

cla(r_ax); plot_looplimit(project.mesh,project.vertices); axis off
hold on; scatter3(project.vertices(:,1),project.vertices(:,2),project.vertices(:,3),'filled')
scatter3(project.vertices(86,1),project.vertices(86,2),project.vertices(86,3),'r','filled')
scatter3(project.vertices(87,1),project.vertices(87,2),project.vertices(87,3),'k','filled')
scatter3(project.vertices(99,1),project.vertices(99,2),project.vertices(99,3),'r','filled')
scatter3(project.vertices(100,1),project.vertices(100,2),project.vertices(100,3),'k','filled')
scatter3(project.vertices(42,1),project.vertices(42,2),project.vertices(42,3),'c','filled')
scatter3(project.vertices(90,1),project.vertices(90,2),project.vertices(90,3),'c','filled')
set(r_ax,'ydir','reverse')

project.images(i).constraints2d = [];
project.images(i).constraints3d = [];
project.images(i).constraintsonsil = [];

while keepGoing
    dcm = datacursormode;
    set(dcm,'Enable','on','DisplayStyle','datatip','SnapToDataVertex','On');
    disp('Select the constraint on the actual dolphin input, then press Return.')
    pause;
    
    cursor = getCursorInfo(dcm);
    project.images(i).constraints2d = vertcat(project.images(i).constraints2d,cursor.Position);
    
    m = input('Does that constraint lie on the silhouette? (y/n)','s');
    if m == 'n'
        project.images(i).constraintsonsil = vertcat(project.images(i).constraintsonsil,[0]);
    else
        project.images(i).constraintsonsil = vertcat(project.images(i).constraintsonsil,[1]);
    end
    
    rotate3d(r_ax,'on');
    dcm = datacursormode;
    set(dcm,'Enable','on','DisplayStyle','datatip','SnapToDataVertex','On');
    disp('Select the corresponding constrain on the dolphin model, then press Return.')
    pause;
    
    rotate3d(r_ax,'off');
    cursor = getCursorInfo(dcm);
    d = sqrt((project.vertices(:,1)-cursor.Position(1)).^2 + (project.vertices(:,2)-cursor.Position(2)).^2 + (project.vertices(:,3)-cursor.Position(3)).^2);
    idx = argmin(d);
    project.images(i).constraints3d = vertcat(project.images(i).constraints3d,idx);

    m = input('Would you like to add another constraint point? (y/n)','s');
    if m ~= 'y'
        keepGoing = false;
    end
end

% update project
if ~onlyConstraints
    project.images(i).rotate = R;
end
project.images(i).translate = [];
project.images(i).transform = [];
project.images(i).scale = [];
nPointsPerSilhouette = floor(length(project.images(i).silhouette)/4);
project.images(i).points = zeros(4,2,nPointsPerSilhouette);

for j=1:4
    idx = [1:4:length(project.images(i).silhouette)-j+1]+j-1;
    idx = idx(1:nPointsPerSilhouette);
    project.images(i).points(j,1,:) = project.images(i).silhouette(idx,1)';
    project.images(i).points(j,2,:) = project.images(i).silhouette(idx,2)';
end

if ~isempty(project.images(i).rotate)
    try
        project = align_project(project,i);
    catch
        disp('Unable to align properly!')
    end
end


%% Check normals

S         = 125;
sil_params        = sil_sample(S, project.images(i).points);
sil_pts           = zeros(S, 2);
sil_normals       = zeros(S, 2);

for s = 1:S
    zerotangent = true;
    while zerotangent
        seg              = floor(sil_params(s));
        sil_pts(s, :) = sil_evalbezier(...
            project.images(i).points(:, :, 1 + seg), ...
            sil_params(s) - seg);
        tan_pts          = 3 * (project.images(i).points(2:end, :, 1 + seg) - ...
            project.images(i).points(1:end - 1, :, 1 + seg));
        tangent          = sil_evalbezier(tan_pts, ...
            sil_params(s) - seg);
        
        zerotangent      = (norm(tangent, 2) == 0);
        sil_params(s)    = sil_params(s) + 1e-4;
    end
    sil_tan = tangent / norm(tangent, 2);
    
    sil_normals(s, :) = [ -sil_tan(2)  sil_tan(1) ];
end

c = 25;
for j=1:length(sil_pts)
    plot(l_ax,[sil_pts(j,1) sil_pts(j,1)+sil_normals(j,1)*c],[sil_pts(j,2) sil_pts(j,2)+sil_normals(j,2)*c],'k','linewidth',1.5)
end

 m = input('Are the normals properly oriented pointing away from the dolphin? (y/n)','s');
 if m ~= 'y'
     project.images(i).normalsLeft = 0;
 else
     project.images(i).normalsLeft = 1;
 end

end
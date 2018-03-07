function label_xyz(ax, s)
% Jason Manley, Sep 2017

if nargin < 1
    ax = gca;  % specify axes to label
end

if nargin < 2
    s = 14;    % specify font size in pixels
end

xlabel(ax,'x','fontsize',s);
ylabel(ax,'y','fontsize',s);
zlabel(ax,'z','fontsize',s);

end
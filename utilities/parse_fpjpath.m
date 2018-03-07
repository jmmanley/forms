function [ points normalsLeft ] = parse_fpjpath(path)

% PARSE_FPJPATH  Parse the silhouette syntax used in .FPJ files
%
%   See also READ_FPJ, READ_PROJECT

assert(path(end) == 'L' || path(end) == 'R');
if path(end) == 'L'
    normalsLeft = true;
else
    normalsLeft = false;
end
path = path(1:end - 1);

assert(path(end) == 'z');
path = path(1:end - 1);

mainpathstart = strfind(path, 'C');
mainpathstart = mainpathstart(1);

startpoint = path(1:mainpathstart - 1);
mainpath = path(mainpathstart:end);

assert(startpoint(1) == 'M');
startpoint = sscanf(startpoint, 'M%f,%f');

points = sscanf(mainpath, 'C%f,%f,%f,%f,%f,%f');
points = reshape(points, 2, 3, []);
points = permute(points, [2 1 3]);
points = [ zeros(1, 2, size(points, 3)) ; points ];

points(1, :, 2:end) = points(end, :, 1:end - 1);
points(1, :, 1) = startpoint;

end


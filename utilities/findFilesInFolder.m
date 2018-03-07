function [files] = findFilesInFolder(folder, fileType, front)
% adapted from Gordon Berman
% Jason Manley, Sep 2017

if nargin==1
    fileType = '.mp4';
end

if nargin < 3 || isempty(front) == 1
    [status,files] = unix(['find ' folder ' -name "*' fileType '"']);
else
    [status,files] = unix(['find ' folder ' -name "' front '*' fileType '"']);
end

if status == 0
    filesOut = regexp(regexp(files,'\t','split'),'\n','split');
    files = filesOut{1}';
    for i=2:length(filesOut)
        files = [files ;filesOut{i}'];
    end
    files = sort(files);
    while isempty(files{1}) == 1
        files = files(2:end);
        if isempty(files)
            break;
        end
    end
    
else
    files = [];
end

end
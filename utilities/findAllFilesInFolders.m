function files = findAllFilesInFolders(folder,fileType,front)
% adapted from Gordon Berman
% Jason Manley, Sep 2017


if nargin==1
    fileType = '.mp4';
end

if nargin < 3
    front = '';
end

[status,files] = unix(['ls ' folder]);
    

if status == 0
    folders = regexp(regexp(files,'\t','split'),'\n','split');
    f = folders{1}';
    for i=2:length(folders)
        f = [f; folders{i}'];
    end
    folders = sort(f);
    while isempty(folders{1}) == 1
        folders = folders(2:end);
    end
    
    files = [];
    for i=1:length(folders)
        currentFiles = findFilesInFolder([folder folders{i}],fileType,front);
        files = [files; currentFiles];
    end
else
    files = [];
end
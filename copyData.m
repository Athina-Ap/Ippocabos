% Source and destination directories
sourceDir = 'Z:\Homes\zutshi01\Recordings\Auditory_Task'; 
destinationDir = 'I:\Data'; 

% List of folders to copy from
foldersToCopy = {'IZ39', 'IZ40', 'IZ43', 'IZ4', 'IZ47', 'IZ48'};

% Loop through each folder
for folderIdx = 1:numel(foldersToCopy)
    currentFolder = foldersToCopy{folderIdx};
    
    sourceFolder = fullfile(sourceDir, currentFolder, 'Final');
    
    % Create the corresponding destination folder
    destinationFolder = fullfile(destinationDir, currentFolder, 'Final');
    mkdir(destinationFolder);
    
    % Copy all files and subfolders except .dat files
    fileList = dir(fullfile(sourceFolder, '*'));
    for fileIdx = 1:numel(fileList)
        if fileList(fileIdx).isdir
            if ~strcmp(fileList(fileIdx).name, '.') && ~strcmp(fileList(fileIdx).name, '..')
                % Recursively copy subfolders
                copyFolderRecursively(fullfile(sourceFolder, fileList(fileIdx).name), fullfile(destinationFolder, fileList(fileIdx).name));
            end
        else
            [~, ~, ext] = fileparts(fileList(fileIdx).name);
            if ~strcmp(ext, '.dat')
                copyfile(fullfile(sourceFolder, fileList(fileIdx).name), destinationFolder);
            end
        end
    end
end

disp('Copy completed.');

function copyFolderRecursively(source, destination)
    mkdir(destination);
    fileList = dir(fullfile(source, '*'));
    for fileIdx = 1:numel(fileList)
        if fileList(fileIdx).isdir
            if ~strcmp(fileList(fileIdx).name, '.') && ~strcmp(fileList(fileIdx).name, '..')
                copyFolderRecursively(fullfile(source, fileList(fileIdx).name), fullfile(destination, fileList(fileIdx).name));
            end
        else
            [~, ~, ext] = fileparts(fileList(fileIdx).name);
            if ~strcmp(ext, '.dat')
                copyfile(fullfile(source, fileList(fileIdx).name), destination);
            end
        end
    end
end

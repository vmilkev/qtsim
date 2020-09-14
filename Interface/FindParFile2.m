function [ patternfiles ] = FindParFile2( changeFolder, pattern )

% Returns a cell array of all files
% with matched pattern

content = dir( changeFolder );

ifile = 0;
idir = 0;

param_file = [];
directories = {};
files = {};

patternfiles = {};

for i = 1:size(content,1)
    if ( content(i).isdir )
        dName = {content(i).name};
        if ( ~( strcmp(dName,'.')||strcmp(dName,'..') ) )
            idir = idir + 1;
            directories(idir,1) = dName;
        end
    else
        ifile = ifile + 1;
        files(ifile,1) = {content(i).name};
    end
end

iFound = 0;

for i =1:size(files,1)
    k = contains( files{i}, pattern );
    if ( k )
        iFound = iFound + 1;
%         name1 = strcat( changeFolder, '/' );
%         name0 = strcat( name1, files{i} );
%         param_file = name0;
        patternfiles{iFound,1} = files{i};
    end
end

if ( isempty(param_file) )
    
    for j = 1:size( directories,1 )
        dirNew = strcat(changeFolder, '/');
        dirNew = strcat(dirNew, directories{j});
        param_file = FindParFile( dirNew, pattern );
        if ( ~isempty(param_file) )
            return;
        end
    end
    
else
    return;
end

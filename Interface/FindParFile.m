function [ param_file ] = FindParFile( changeFolder, pattern )

content = dir( changeFolder );

ifile = 0;
idir = 0;

param_file = [];
directories = {};
files = {};

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

for i =1:size(files,1)
    k = contains( files{i}, pattern );
    if ( k )
        name1 = strcat( changeFolder, '/' );
        name0 = strcat( name1, files{i} );
        param_file = name0;
        return
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


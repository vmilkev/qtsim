function [ res ] = MakeFolder( dirName )

status = [];
msg = [];

res = 1;

if ~exist(dirName, 'dir')
    [ status, msg ] = mkdir(dirName);
end

if (~status)
    runtime_log( 0, 'ERROR: cannot create a folder to store simulation results!' );
    runtime_log( 0, msg );
    res = 0;
end

if ( ~isempty(msg) )
    runtime_log( 0, msg );
end

end
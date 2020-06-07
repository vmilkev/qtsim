function runtime_log( new_sim, message, txt_value )

if ~exist('txt_value','var')
    txt_value = [];
end

% res_folder = strcat( 'traitsim_results-', date );
% logFileName = strcat( res_folder, '/runtime.log' );
logFileName = strcat( 'tsim_runtime.log' );

% if ( ~MakeFolder(res_folder) )
%     return;
% end

if (new_sim)
    fileID = fopen(logFileName,'w');
else
    fileID = fopen(logFileName,'a');
end

if ( isempty(txt_value) )
    fprintf(fileID,'%12s %2s %12s\n',datetime('now'), '| ', message);
else
    fprintf(fileID,'%12s %2s %12s %12s\n',datetime('now'), '| ', message, txt_value);
end

fclose(fileID);

end
function [ P ] = read_products( prod_files )

%     prod_files is a cell array consisting
%     of products files names

P = read_prod_file( prod_files{1,1} );

if ( size(prod_files,1) > 1 )
    for i = 2:size(prod_files,1)
        tP = read_prod_file( prod_files{i,1} );
        P = [ P; tP(3:end,:) ];
    end
end


end


% LOCAL FUNCTION

function [ P ] = read_prod_file( prodFileName )

fileID = fopen(prodFileName,'r');

if ( fileID == -1 )
    runtime_log( 0, 'ERROR: cannot open product file!' );
    return;
end

line = 0;
i = 1;

while ( line ~= -1 )
    line = fgetl(fileID);
    
    if ( line == -1 )
        break;
    end
    
    C = textscan(line,'%f');
    P(:,i) = C{1,1};

    i = i + 1;
end

fclose(fileID);

end
function [values, counts] = degree_distr(...
                                        d,...
                                        fname,...
                                        printInFile...
                                        )

% Function which calculates connectivity (graph degree) distribution.

% d - network connectivity (graph degree);
% fname - parameter of type string representing a file name where the
%         vectors values and counts will be printed if printInFile = TRUE;
% printInFile - logical parameter, if printInFile = TRUE data will be printed to a specified file;
% values - vector of gene's connectivity values (graph node's degrees);
% counts - vector of connectivity frequencies multiplied by a total number of genes.

valNum = numel(d);
vMax = max(d);
degree = 1;
j = 1;
freeFlag = true;
while (degree <= vMax)
    for i = 1:valNum
        if (d(i,1) == degree)
            if (freeFlag)
                counts(j,1) = 0;
            end
            counts(j,1) = counts(j,1) + 1;
            freeFlag = false;
        end
    end
    if (~freeFlag)
        values(j,1) = degree;
        j = j + 1;
        freeFlag = true;
    end
    degree = degree + 1;
end

A(:,1) = values(:,1);
A(:,2) = counts(:,1);

if (printInFile)
    fileID = fopen(fname,'w');
    fprintf(fileID,'%10s %10s\n','degree','counts');
    fprintf(fileID,'%10d %10d\n', A');
    fclose(fileID);
end

end
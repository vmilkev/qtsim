function [ SNP, indID ] = ImportSNPtoMatrix( filename )

SNP = [];
indID = [];

T = readtable( filename );
temp_snp = table2array(T(:,1:end));

snpSize1 = size(temp_snp);
indID = zeros(1,snpSize1(1,1));
indID(1,:) = temp_snp(:,1);

SNP = zeros( snpSize1(1,1), snpSize1(1,2)-2 );
SNP(:,:) = temp_snp(:,3:end);

snpSize = size(SNP);

% find missing data (values other than 0,1,2) and sabstitute to
% alleles with highest frequencies values

for i = 1:snpSize(1,2)
    
    p0 = 0;
    p1 = 0;
    p2 = 0;
    counts = 0;
    for j = 1:snpSize(1,1)
        all = SNP(j,i);
        switch all
            case 0
                p0 = p0 + 1;
                counts = counts + 1;
            case 1
                p1 = p1 + 1;
                counts = counts + 1;
            case 2
                p2 = p2 + 1;
                counts = counts + 1;
        end
    end
    p0 = p0/counts;
    p1 = p1/counts;
    p2 = p2/counts;
    
    pTmp = 0;
    if (p0 >= p1)
        maxAll = 0;
        pTmp = p0;
    else
        maxAll = 1;
        pTmp = p1;
    end
    if (p2 >= pTmp)
        maxAll = 2;
    end
    
    for j = 1:snpSize(1,1)
        all = SNP(j,i);
        if ( all > 2 )
            SNP(j,i) = maxAll;
        end
    end
    
end

end
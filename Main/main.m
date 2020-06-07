function main( setupfile )

warning('off', 'MATLAB:datetime:NonstandardSystemTimeZone');

pool_1 = gcp;

if ~exist('setupfile','var')
    tsimfile = [];
else
    tsimfile = setupfile;
end

RunTraitSimulator( tsimfile );

delete(pool_1);

end

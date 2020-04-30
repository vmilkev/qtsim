function main()

warning('off', 'MATLAB:datetime:NonstandardSystemTimeZone');

pool_1 = gcp;

RunTraitSimulator();

delete(pool_1);

end

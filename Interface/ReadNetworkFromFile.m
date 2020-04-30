function [ N ] = ReadNetworkFromFile( fileName )

%N = dlmread( fileName, ' ');
N = full(spconvert(load( fileName )));

end
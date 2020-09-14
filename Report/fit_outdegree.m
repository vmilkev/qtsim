function [a,b,rsquare,xPred,yPred] = fit_outdegree(...
                                                  degree,...
                                                  counts...
                                                  )

% Approximation of out-Connectivity distribution data (degree & counts) by
% power-law fit.
% [a,b,rsquare,xPred,yPred] - parameters of yPred =  a*xPred.^b;

[xData, yData] = prepareCurveData( degree, counts );

ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';

% % Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

a = fitresult.a;
b = fitresult.b;
rsquare = gof.rsquare;
xPred = linspace(1,max(degree),100);
yPred = a*xPred.^b;

end
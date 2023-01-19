function [yMean, yCI95] = mean_confidence_interval(y)
%Input: data (number of samples x independent variable length)
%
%Output: mean and conficence intervals

y = y.';
N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
yMean = mean(y);                                    % Mean Of All Experiments At Each Value Of the independent variable
ySEM = std(y)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of the independent variable
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of the independent variable

yMean = yMean.';
yCI95 = yCI95.';
end
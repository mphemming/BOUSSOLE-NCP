%% This file shows a couple examples of the error MC propagation 
% Carsten Robens and Stefan Brakhane
% Robens@iap.uni-bonn.de Brakhane@iap.uni-bonn.de
% Created: 11.05.2016
% Last Modified: 25.05.2016

%% Simple example to show that gaussian Error propagation holds true for linear functions
A = generateMCparameters('gaussian',[2,0.2]);
B = generateMCparameters('gaussian',[0.5,0.2]);
paramMatrix = [A;B];
funToProp = @(x) x(1)+x(2);
[funValue,funCI,funSamples] = propagateErrorWithMC(funToProp, paramMatrix);

fprintf('Traditional Gaussian Error Propagartion: %.2f [%.2f - %.2f]\n',(2+0.5),(2+0.5)-sqrt(0.2^2+0.2^2),(2+0.5)+sqrt(0.2^2+0.2^2));

%% Simple example to show that divisions can lead to assymetric distrubution
A = generateMCparameters('gaussian',[2,0.2]);
B = generateMCparameters('gaussian',[0.5,0.2]);
paramMatrix = [A;B];
funToProp = @(x) x(1)./x(2);
[funValue,funCI,funSamples] = propagateErrorWithMC(funToProp, paramMatrix);

fprintf('Traditional Gaussian Error Propagartion: %.2f [%.2f - %.2f]\n',(2/0.5),(2/0.5)-sqrt((1/0.5)^2*0.2^2+(2/(0.5^2))*0.2^2),(2/0.5)+sqrt((1/0.5)^2*0.2^2+(2/(0.5^2))*0.2^2));

%% since this example shows a clear assymetrie we should also consider what we define as central function value
[funValue,funCI,funSamples] = propagateErrorWithMC(funToProp, paramMatrix,'method','mean');

%% since this example shows a clear assymetrie we should also consider what we define as central function value
[funValue,funCI,funSamples] = propagateErrorWithMC(funToProp, paramMatrix,'method','maximum');

%% combining different initial error distributions (gaussian and binomial)
A = generateMCparameters('gaussian',[1,0.02],'plot',true);
B = generateMCparameters('binomial',[20,18],'plot',true);
paramMatrix = [A;B];
funToProp = @(x) x(1).*x(2);
[funValue,funCI,funSamples] = propagateErrorWithMC(funToProp, paramMatrix);


%% combining different initial error distributions (gaussian and bootstrapping Mean)
A = generateMCparameters('gaussian',[5,0.2],'plot',true);
B = generateMCparameters('bootstrapMean',[5.9,6.8,6.1,5.8,5.6,5.7,6.1,6.6,6.2,6.1,5.5],'plot',true);
paramMatrix = [A;B];
funToProp = @(x) x(1)*x(2);
[funValue,funCI,funSamples] = propagateErrorWithMC(funToProp, paramMatrix);


%% combining different initial error distributions (gaussian and bootstrapping Distribution)
A = generateMCparameters('gaussian',[5,0.2],'plot',true);
B = generateMCparameters('bootstrapDistribution',[5,7,6,5,5,5,6,7,6,6,5],'plot',true);
paramMatrix = [A;B];
funToProp = @(x) x(1)*x(2);
[funValue,funCI,funSamples] = propagateErrorWithMC(funToProp, paramMatrix);


%% you can also use the return of the error propagation for another function
A = generateMCparameters('gaussian',[2.0,0.2]);
B = generateMCparameters('gaussian',[1,0.2]);

paramMatrix = [A;B];
funToProp = @(x) (x(1)+x(2))./abs(x(1)-x(2));
[funValue,funCI,funSamples] = propagateErrorWithMC(funToProp, paramMatrix);

C = generateMCparameters('gaussian',[3,0.2]);
funToProp = @(x) x(1)*x(2);
paramMatrix = [funSamples;C];
[funValue,funCI,funSamples] = propagateErrorWithMC(funToProp, paramMatrix);


%% Calculating the average n from a sideband ratio measurement
A = generateMCparameters('gaussian',[0.033,0.021]);
B = generateMCparameters('gaussian',[0.542,0.052]);
paramMatrix = [A;B];
funToProp = @(x) (x(1)./x(2))./(1-(x(1)./x(2)));
[funValue,funCI,funSamples] = propagateErrorWithMC(funToProp, paramMatrix);
fprintf('Average Vibrational Excitation nbar:  %.3f [%.3f - %.3f]\n',funValue,funCI(1),funCI(2));

%% and the resulting groundstate population
paramMatrix = [funSamples];
funToProp = @(x) 1./(1+x(1));
[funValue,funCI,~] = propagateErrorWithMC(funToProp, paramMatrix);
fprintf('Groundstate Occupation n0:  %.3f [%.3f - %.3f]\n',funValue,funCI(1),funCI(2));
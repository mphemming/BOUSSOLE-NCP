function [funValue, funCI, mcSamples] = propagateErrorWithMC(funOfInterest, params, varargin)
%propagateErrorWithMC 
%
% Inputs:
%   funOfInterest  function that should be evaluated
%   params         matrix of column vectors, each row represents sampled parameters
%   CIthreshold    (optional) confidence interval threshold, default 0.68
%   plot           (optional) plot final distribution
%   method         (optional) method to determine funValue (median (default), mean, maximum)
%  
% Outputs:
%   funValue   calculated function value
%   funCI      confidence interval for function value
%   mcSamples  samples of the mc function distribution

% Carsten Robens and Stefan Brakhane
% Robens@iap.uni-bonn.de Brakhane@iap.uni-bonn.de
% Created: 11.05.2016
% Last Modified: 17.05.2016

    p = inputParser();
    p.addParameter('CIthreshold', 0.68, ...
                   @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
    p.addParameter('plot', true, ...
                   @islogical);
    p.addParameter('method', 'median', ...
                   @(x) any(strcmp(x,{'median', 'mean', 'maximum'}))); 

    p.parse(varargin{:});

    CIthreshold = p.Results.CIthreshold;
    MonteCarlo_length = round(size(params,2)/100)*100; % need to be a multiple of 100
    funVals = zeros(1,MonteCarlo_length);
    
    fprintf('Evaluating function values: %3.0f%%\n',0);
    
    for idx = 1:100
        for kdx = 1:ceil(MonteCarlo_length/100)
            funValIdx = (idx-1)*ceil(MonteCarlo_length/100)+kdx;
            funVals(funValIdx) = funOfInterest(params(:,funValIdx));
        end
        fprintf('\b\b\b\b\b%3.0f%%\n',idx);
    end
    

    %% calculate interval
    funError = std(funVals);
    funMean = mean(funVals);

    funBins = linspace(funMean-funError, funMean+funError, 20).';
    vals_bins = histc(funVals,funBins);
    
    delLeft  = sum(funVals<funBins(1));
    delRight = sum(funVals>funBins(end));
    
    cmax = MonteCarlo_length/1000;
    iter = 0;
    numIterBinAdjust = 20;
    
    fprintf('Adjusting bin width:        %3.0f%%\n',0);
    while iter<numIterBinAdjust
        if delLeft>cmax
            funBins = linspace(funBins(1)-funError, funBins(end), numel(funBins) + 10);
%             fprintf('extend left ');
        end
        
        if delRight>cmax
            funBins = linspace(funBins(1), funBins(end)+funError, numel(funBins) + 10);
        end
        
        if delLeft<cmax && delRight<cmax
            break;
        end
        
        vals_bins = histc(funVals,funBins);
        delLeft  = sum(funVals<funBins(1));
        delRight = sum(funVals>funBins(end));
        
        iter = iter + 1;
        fprintf('\b\b\b\b\b%3.0f%%\n',iter/numIterBinAdjust*100);
    end
    fprintf('\b\b\b\b\b%3.0f%%\n',100);
    
    %% determine bin width
    iter = 0;
    maxNumOfBinAdjust=20*20/length(funBins);
    fprintf('Adjusting number of bins:   %3.0f%%\n',0);
    while iter<maxNumOfBinAdjust && max(diff(vals_bins/sum(vals_bins)))>0.01
        funBins = linspace(funBins(1), funBins(end), numel(funBins)*2);
        vals_bins = histc(funVals, funBins);
        iter = iter + 1;
        fprintf('\b\b\b\b\b%3.0f%%\n',iter/maxNumOfBinAdjust*100);
    end
    fprintf('\b\b\b\b\b%3.0f%%\n',100);
    
    %% normalize
    vals_bins = vals_bins(1:end-1)/sum(vals_bins);
    vals_bins_centers = 0.5 * (funBins(2:end) + funBins(1:end-1));
    
    mcSamples = funVals;

    
    %% confidence interval
    cumSumLowerPart = cumtrapz(vals_bins_centers, vals_bins);
    cumSumLowerPart = cumSumLowerPart/max(cumSumLowerPart);
    
    mask = [false, diff(cumSumLowerPart)==0];
    cumSumMask = cumSumLowerPart<0.99 & cumSumLowerPart>0.01 & ~mask;

    lowerCI_outwardsIntegration = interp1(cumSumLowerPart(cumSumMask),vals_bins_centers(cumSumMask),(1-CIthreshold)/2);
    upperCI_outwardsIntegration = interp1(1-(cumSumLowerPart(cumSumMask)),(vals_bins_centers(cumSumMask)),(1-CIthreshold)/2);
    
    funCI = [lowerCI_outwardsIntegration, upperCI_outwardsIntegration];
    
    %% function value
    switch(p.Results.method)
        case 'median'       % determined by median of function values
            funValue = median(funVals);
        case 'mean'  % value as expectation value
            funValue = mean(funVals);
        case 'maximum'  % value as maximum of probability distribution
            [~, id] = max(vals_bins);
            sel = max([1 id-4]) : min([numel(vals_bins), id+4]);
            pf = polyfit(vals_bins_centers(sel), vals_bins(sel), 2);
            fu = @(x) polyval(pf, x);

            [~, id] = max(vals_bins);
            funValue = fminsearch(@(x) -fu(x), vals_bins_centers(id));
    end
    
    
    %% plotting
    if p.Results.plot
        figure;
        idxP = find(vals_bins*MonteCarlo_length>20, 1, 'first') : find(vals_bins*MonteCarlo_length>20, 1, 'last');
        
        CImask = vals_bins_centers(idxP)>funCI(1) & vals_bins_centers(idxP)<funCI(2);
        bar(vals_bins_centers(idxP(CImask)), vals_bins(idxP(CImask)), 0.5,'FaceColor',[49,163,84]/255,'EdgeColor',[0,109,44]/255,'LineWidth',1);
        hold on;
        bar(vals_bins_centers(idxP(~CImask)), vals_bins(idxP(~CImask)), 0.5,'FaceColor',[49,130,189]/255,'EdgeColor',[8,81,156]/255,'LineWidth',1);
        xlabel('Function Values', 'FontSize', 16);
        ylabel('Montecarlo Stastistical Distribution', 'FontSize', 16);

        gauss_fit_funct = @(p,x) p(1)*exp(-0.5*(x-p(2)).^2/p(3).^2);
        pFit = nlinfit(vals_bins_centers, vals_bins, gauss_fit_funct, [0.2, funValue, (funCI(2)-funCI(1))/2]);
        
        plot(vals_bins_centers(idxP), gauss_fit_funct(pFit,vals_bins_centers(idxP)), 'Color',[255,51,51]/255, 'LineWidth',3);
        
        if strcmp(p.Results.method,'maximum')
            plot(vals_bins_centers(sel), fu(vals_bins_centers(sel)), 'Color','g', 'LineWidth',3);
        end
        y1=get(gca,'ylim');
        plot([funValue funValue],y1,'--','color','k','LineWidth',1.5)
        hold off
        drawnow;

        fprintf('Monte Carlo:  %.2f [%.2f - %.2f]\n', funValue, funCI(1), funCI(2));
        fprintf('Gaussian Fit: %.2f [%.2f - %.2f]\n', pFit(2), pFit(2)-abs(pFit(3)), pFit(2)+abs(pFit(3)));
    end
end
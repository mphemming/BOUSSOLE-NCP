function [smooth_struct] = LOESS_window(time, var, pm_day)

%% remove NaNs

check_nans = isnan(time) & isnan(var);
time(check_nans) = [];
var(check_nans) = [];

%% identify days availabe in time
[yrs,mns,dys,~,~,~] = datevec(time);
avail_days = unique(datenum(yrs,mns,dys));
interval = nanmean(diff(time));
span = (pm_day*2)/interval;
span = span/numel(var);

%% Smooth the data points
% define variables
x = 1:1:length(time);
y = var;
n = length(x);
r = x(end) - x(1);
hlims = [span,x(1);... 
    (span)/2,x(1)+r*span/2;...
    (span)/2,x(1)+r*(1-span/2);...
    span,x(end)];  
% Find the LOESS fit to the data
smoothed = zeros(n,1); % pre-allocate space
for i = 1:n
    
    % define the tricube weight function
    h = interp1(hlims(:,2),hlims(:,1),x(i));
    % w = (1-abs((x./max(x)-x(i)./max(x))/h).^3).^3;
    w = ones(size(time)); % no weighting
    
    % data points outwith the defined span can be ignored (for speed)
    w_idx = w>0;
    w_ = w(w_idx);
    x_ = x(w_idx);
    y_ = y(w_idx);
    
    % Calculate the weighted coefficients
    XX =   [nansum(w_.*x_.^0), nansum(w_.*x_.^1), nansum(w_.*x_.^2);...
            nansum(w_.*x_.^1), nansum(w_.*x_.^2), nansum(w_.*x_.^3);...
            nansum(w_.*x_.^2), nansum(w_.*x_.^3), nansum(w_.*x_.^4)];
    
    YY =   [nansum(w_.*y_.*(x_.^0));...
            nansum(w_.*y_.*(x_.^1));...
            nansum(w_.*y_.*(x_.^2))];
    
    CC = XX\YY;
    
    % calculate the fitted data point
    smoothed(i) = CC(1) + CC(2)*x(i) + CC(3)*x(i).^2;
    
    % add var in each span
    smooth_struct.calc(i).times = time(x_);
    smooth_struct.calc(i).O2 = y_;
    smooth_struct.calc(i).offset = CC(1);
    smooth_struct.calc(i).slope_1= CC(2);
    smooth_struct.calc(i).slope_2= CC(3);
    [~,mns,dys,~,~,~] = datevec(smooth_struct.calc(i).times);
    smooth_struct.calc(i).unique_days = unique(dys); 
    smooth_struct.calc(i).unique_days_n = numel(unique(dys)); 
end

%% save information
smooth_struct.smoothed = smoothed;
smooth_struct.time = time;
smooth_struct.var = var;
smooth_struct.span = span;
smooth_struct.pm_day = pm_day

end
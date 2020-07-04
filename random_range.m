function r = random_range(n_min, n_max,n_samples)

r = (n_max-n_min).*rand(n_samples,1) + n_min;
r(1:numel(r)/2) = r(1:numel(r)/2)*-1;

end
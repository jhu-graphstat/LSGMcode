function profile_likelihood = profileLikelihood(d)
n = length(d);

profile_likelihood = zeros(1,n);

for q = 1:n
    sample_1 = d(1:q);
    sample_2 = d((q+1):end);

    mu_1_hat = mean(sample_1);
    mu_2_hat = mean(sample_2);
    sigma_hat = ((q-1)*var(sample_1) + (n-q-1)*var(sample_2))/(n-2);

    pl_sample_1 = -(sample_1 - mu_1_hat).^2/(2*sigma_hat) - 0.5*log(2*pi*sigma_hat);
    pl_sample_2 = -(sample_2 - mu_2_hat).^2/(2*sigma_hat) - 0.5*log(2*pi*sigma_hat);

    profile_likelihood(q) = sum(pl_sample_1) + sum(pl_sample_2);
end
end
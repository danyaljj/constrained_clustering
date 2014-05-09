function out = gaussian_log_pdf(x, mu, sigma)

% check determinant 
det1 = det(sigma); 
if( det1 < 0 ) 
    error('Oh no! determinant is singular!!')
end 

inv_sigma = inv(sigma);
log_norm_const = -0.5 * log(det1) - size(x,2) / 2.0 * log(2*pi);
log_result = -0.5*( (x - mu) * inv_sigma * (x - mu)' );
out = log_norm_const + log_result;

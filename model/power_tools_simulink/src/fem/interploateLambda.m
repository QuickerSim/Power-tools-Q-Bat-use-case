function lambda_interp = interploateLambda(lambda, u, flag_is_cylindrical)
coder.allowpcode('plain')

field_names = fieldnames(lambda);
% count number of lambdas occurence
n_lambdas = 0;
for i = 1 : size(field_names, 1)
    if count(field_names{i}, 'lambda')
        n_lambdas = n_lambdas + 1;
    end
end

if(n_lambdas == 3)
    if flag_is_cylindrical
        
        lambda_interp = zeros(3, length(u));
        
        % phi component
        lambda_interp(1, :) = nearestInterp1(lambda.T, lambda.lambda_x, u);
        % r component
        lambda_interp(2, :) = nearestInterp1(lambda.T, lambda.lambda_y, u);
        % z component
        lambda_interp(3, :) = nearestInterp1(lambda.T, lambda.lambda_z, u);
        
    else
        
        lambda_interp = zeros(3, 3, length(u));
        
        % x component?
        lambda_interp(1, 1, :) = nearestInterp1(lambda.T, lambda.lambda_x, u);
        % y component?
        lambda_interp(2, 2, :) = nearestInterp1(lambda.T, lambda.lambda_y, u);
        % z component?
        lambda_interp(3, 3, :) = nearestInterp1(lambda.T, lambda.lambda_z, u);
    end
    
elseif(n_lambdas == 1)
    
    lambda_interp = nearestInterp1(lambda.T, lambda.lambda, u);
    
else
    
    if size(lambda.lambda, 1) == 1
        lambda_interp = lambda * ones(length(u), 1);
    else
        lambda_interp = repmat(lambda.lambda, 1, 1, length(u));
    end
    
end

end
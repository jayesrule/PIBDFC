%% function to create transition matrix
function trans_mat = create_trans_mat(xi_est,  rho_est, x_level)
    [S,~] = size(rho_est);
    
    oomph_add = zeros(size(xi_est));
    
    for s = 1:S
        oomph_add(:,s) = rho_est(s,:)*x_level;
    end
    
    trans_mat  = xi_est+oomph_add;
    trans_mat = exp(trans_mat)./repmat(sum(exp(trans_mat),2),[1,S]);
    
    

end
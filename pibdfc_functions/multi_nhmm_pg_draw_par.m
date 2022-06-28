function [sqrl_all, Z, rho_all, eta] = multi_nhmm_pg_draw_par(xDat, sqrl_old, Z_old, Z_prior, rho_old, eta_old, eta_prior, s_assign, a0, b0)

%% This function is a Gibbs sampling step to sample new values for the log linear coefficients
%% in the Holsclaw et al. paper for sampling from heterogenous HMMs.

%% Inputs:
% Xi_old: The values of the Xi parameter for the previous iteration
% Z: The arrangement of the state transitions in binary form for the current iteration
% Exo: The exogenous variables incorporated in the log linear regression, X in the Holsclaw paper
% a: prior mean of the Xi, can be scalar or assign unique prior value to each term
% b: prior variance of the Xi, can be scalar or assign unique prior value to each term

%% Outputs:
% Xi_new: The sampled values of the Xi parameter for the current iteration


%%
[T,K,N] = size(xDat);
sqrl_all = sqrl_old;
Z = Z_old;
rho_all = rho_old;
eta = eta_old;
S = size(Z,1);

s_lag = [s_assign(2:end,:);zeros(1,N)];

%Update Subject level sqrl for each subject
parfor i = 1:N
    
    sqrl = sqrl_all(:,:,i);
    rho = rho_all(:,:,i);
    %Compute transition count matrices
    trans_counts = zeros(S, S);
    for r = 1:S
        for s = 1:S
            trans_counts(r,s) = sum((s_assign(:,i)==r) & (s_lag(:,i)==s));
        end
    end
    
    trans_from = sum(trans_counts,1);
    
    %updating sqrl
    for j = 1:S
        for k = 2:S
            omega = zeros(T,1);
            c_jki = zeros(T,1);
            for t = 1:T
                c_jki(t) = -squeeze(xDat(t,:,i))*squeeze(rho(k,:))'+log(sum(exp(sqrl(j,setdiff(1:S,k))).*exp(squeeze(xDat(t,:,i))*squeeze(rho(setdiff(1:S,k),:))')));
                pg_scale = sqrl(j,k) - c_jki(t);
                omega(t) = pgdraw(pg_scale);
            end
            m_omega = Z(j,k)/a0+omega'*c_jki+trans_counts(j,k)-trans_from(j)/2;
            m_omega = m_omega/(1/a0+sum(omega));
            v_omega = 1/sqrt(sum(omega)+1/a0);
            %sqrl(j,k,i) = max(-10,min(10,normrnd(m_omega, v_omega)));
            sqrl(j,k) = normrnd(m_omega, v_omega);
        end
    end
    
    %updating rho
    if K > 1
        for j = 2:S
            psi = 1*(1*s_assign(:,i) == j);
            psi = [psi(2:end);0];
            for k = 1:K
                omega = zeros(T,1);
                c_jki = zeros(T,1);
                for t = 1:T
                    s = s_assign(t,i);
                    cov_set = setdiff(1:K,k);
                    state_set = setdiff(1:S,j);
                    c_jki(t) = -sqrl(s,j)-squeeze(xDat(t,cov_set,i))*squeeze(rho(j,cov_set))+log(sum(exp(sqrl(s,state_set)).*exp(xDat(t,:,i)*squeeze(rho(state_set,:))')));
                    pg_scale = xDat(t,k,i)*rho(j,k) - c_jki(t);
                    omega(t) = pgdraw(pg_scale);
                end
                m_omega = eta(j,k)/a0+sum((psi-0.5+omega.*c_jki).*xDat(:,k,i));
                m_omega = m_omega/(1/a0+sum(omega.*(xDat(:,k,i).^2)));
                v_omega = 1/sqrt(1/a0+sum(omega.*(xDat(:,k,i).^2)));
                rho(j,k) = normrnd(m_omega, v_omega);
            end
        end
        sqrl_all(:,:,i) = sqrl;
        rho_all(:,:,i) = rho;
        
        
    else
        for j = 2:S
            k = 1;
            psi = 1*(1*s_assign(:,i) == j);
            psi = [psi(2:end);0];
            omega = zeros(T,1);
            c_jki = zeros(T,1);
            for t = 1:T
                s = s_assign(t,i);
                state_set = setdiff(1:S,j);
                c_jki(t) = -sqrl(s,j)+log(sum(exp(sqrl(s,state_set)).*exp(xDat(t,:,i)*squeeze(rho(state_set,:))')));
                pg_scale = xDat(t,k,i)*rho(j,k) - c_jki(t);
                omega(t) = pgdraw(pg_scale);
            end
            m_omega = eta(j,k)/a0+sum((psi-0.5+omega.*c_jki).*xDat(:,k,i));
            m_omega = m_omega/(1/a0+sum(omega.*(xDat(:,k,i).^2)));
            v_omega = 1/sqrt(1/a0+sum(omega.*(xDat(:,k,i).^2)));
            rho(j,k) = normrnd(m_omega, v_omega);
        end
    end
    
    sqrl_all(:,:,i) = sqrl;
    rho_all(:,:,i) = rho;
    
end





%update global params
sqrl_bar = mean(sqrl_all,3);
Z = (b0*N*sqrl_bar+a0*Z_prior)/(a0+b0*N) + sqrt(a0*b0/(a0+N*b0))*normrnd(0,1,[S,S]);
Z(:,1) = 0;

rho_bar = mean(rho_all,3);
eta = (b0*N*rho_bar+a0*eta_prior)/(a0+b0*N) + sqrt(a0*b0/(a0+N*b0))*normrnd(0,1,[S,K]);
eta(1,:) = 0;


%%


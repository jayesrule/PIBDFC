function [map_states] = mapstate_fb_draw(post_draws, Y, xDat)
%%
st_post = post_draws.states_post;
Omega_est = mean(post_draws.Omega_post,4);
sqrl = mean(post_draws.sqrl_post,4);
rho = mean(post_draws.rho_post,4);

N = size(st_post,1);
subs = size(Y,3);
S = size(Omega_est,3);
R = size(Y,2);
map_states = zeros(N,subs);

%%
for sub = 1:subs
    h = zeros(N,1);
    L = zeros(N,S);
    r = st_post(1,sub,1);
    L(1,r) = mvnloglike(squeeze(Y(1,:,sub)), repelem(0,R),Omega_est(:,:,r));
    for n = 2:N
        for s = 1:S
            temp = zeros(S,1);
            for r = 1:S
                q = sqrl(r,s,sub)+squeeze(xDat(n,:,sub))*squeeze(rho(s,:,sub))'-log(sum(exp(sqrl(r,:,sub)+squeeze(xDat(n,:,sub))*squeeze(rho(:,:,sub))')));
                temp(r) = L(n-1,r)+q;
            end
            L(n,s) = max(temp) + mvnloglike(squeeze(Y(n,:,sub)), repelem(0,R),Omega_est(:,:,s));
        end
    end
    
    [~,h(N)] = max(L(N,:));
    
    for n = flip(1:(N-1))
        q = sqrl(:,h(n+1),sub)+squeeze(xDat(n,:,sub))*squeeze(rho(h(n+1),:,sub))'-log(sum(exp(sqrl(:,:,sub)+squeeze(xDat(n,:,sub))*squeeze(rho(:,:,sub))')))';
        [~,h(n)] = max(L(n,:)+q');
    end
    
    map_states(:,sub) = h;

end

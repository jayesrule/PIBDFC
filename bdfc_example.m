%This script will generate data and then fit pibdfc as an example.
%%
addpath('pibdfc_functions');
addpath('bdfc_functions');
rng(34);

%% Generating the data
run('create_gaussian_simset.m');

%% Concatonating dataset
y_dat = zeros(N*Subs,R);
true_states_long = reshape(true_states(:,1:Subs),N*Subs,1);

for i = 1:Subs
    y_dat((1+(i-1)*N):(i*N),:) = Y(:,:,i);
end
xDat = zeros(N*Subs,1);

%% Setting Up for the MCMC
nsims = 1000; %Number of simulations total
nburn = 1000; %Number of burn-in sims
S = 3;
fprintf('Fitting bdfc: \n')
%Perform MCMC
post_draws = Call_Function(y_dat,xDat,S,nburn,nsims,34);

%% Relabeling
map_states_long = mapstates(post_draws.states');
[map_states_long, map] = relabel_results2(map_states_long', true_states_long);
map_states = reshape(map_states_long', N,Subs);

%% Creating summaries
est_tot = bdfc_fdr_thresh(post_draws,0.05, false);
est_tot = est_tot(:,:,map);
parcor = prec2parcor(mean(post_draws.C_save,4));
parcor = parcor(:,:,map);
parcor(est_tot==0) = 0;

%% Plot partial correlation matrices

%Plotting the true partial correlation on top
for i=1:3
    subplot(2,S,i)
    h = heatmap(prec2parcor(Omegas(:,:,i)));
    h.ColorLimits = [-1 1];
    h.Colormap = bluewhitered;
    title(sprintf('Truth: State %d',i));
    if i<3
        h.ColorbarVisible = 'off';
    else
        annotation('textarrow',[1,1],[0.75,0.75],'string','Partial Correlation', ...
            'HeadStyle','none','LineStyle','none','HorizontalAlignment','center','TextRotation',90);
    end
end



%PLotting the estimated partial correlation on the bottom
for i=1:3
    subplot(2,S,i+S)
    h = heatmap(parcor(:,:,i));
    h.ColorLimits = [-1 1];
    h.Colormap = bluewhitered;
    title(sprintf('PIBDFC: State %d',i));
    if i<3
        h.ColorbarVisible = 'off';
    else
        annotation('textarrow',[1,1],[0.25,0.25],'string','Partial Correlation', ...
            'HeadStyle','none','LineStyle','none','HorizontalAlignment','center','TextRotation',90);
    end
end

%% State Sequence Plot

%RGB for each state
Viterbi_colors = [243 118 0; 190 248 232;  132 195 88;...
    231 119 153; 182 108 108; 145 132 233]/255;
Viterbi_choose = Viterbi_colors(1:S,:); %Choosing the first S from Viterbi_colors



figure()
axes('Position', [0.08, 0.54, 0.83 ,0.38]);
customSubLabel = string(1:Subs);
customSubLabel(1:end) = '';

%Cluster similar state trajectories together
tree = linkage(true_states', 'average');
[H,~,outperm] = dendrogram(tree,'Orientation','left', 'Labels',cellstr(customSubLabel));
outperm = flip(outperm);
set(gca,'XColor','none');
set(gca,'YColor','none');

%Create image
h = image(true_states(:,outperm)');
colormap(Viterbi_choose);
hold on;
for K = 1 : 3; hidden_h(K) = surf(uint8(K-[1 1;1 1]), 'edgecolor', 'none'); end
hold off
uistack(hidden_h, 'bottom');
lgd = legend(hidden_h, {'1', '2', '3'});
title(lgd,'Conn State')
xline(301, 'k--','LineWidth', 1, ...
    'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','bottom',...
    'Alpha',1, 'FontSize', 10)

lgd.String = {'1' '2' '3'};
lgd.Position = [0.92,0.5,0.05,0.05];
title('True State Transition Paths')
xticks([1,(1:21)*50])
ylabel('Subject Index');





axes('Position', [0.08, 0.08, 0.83 ,0.38]);
h = image(map_states(:,outperm)');
colormap(Viterbi_choose);
xline(301, 'k--','LineWidth', 1, ...
    'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','bottom',...
    'Alpha',1, 'FontSize', 10)

title('PIBDFC Estimated State Transition Path by Subject')
xticks([1,(1:21)*50])
xlabel('TR');
ylabel('Subject Index');


%% Change Point Plot for subject: sub
sub = 1; %Subject to inspect
inds = (1+(sub-1)*N):(sub*N);
inds = inds(2:end);

c_point_truth = zeros(1, N);
c_point_est = c_point_truth;
c_point_truth(2:end) = 1*(true_states(2:N,sub) ~= true_states(1:(N-1),sub))';
c_point_est(2:end) = mean(post_draws.states(inds,:) ~= post_draws.states(inds-1,:),2)';
figure()
plot(1:N, c_point_truth(sub,:), 'k-', 'LineWidth',1.5);
ylim([0,1]);
title('Change Point Estimation for Sim Subject 1')
xlabel('TR')
ylabel('P(s_t \neq s_{t+1}|Y^1)')
hold on
plot(1:N, c_point_est(sub,:),'r-');
yline(0.95,'k--')
lgd = legend('True Change Point', 'Est of Changepoint');
ldg.Location = 'northeastoutside';


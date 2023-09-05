clear all; clc; close all; addpath('functions\'); rng('default')
%% USER INPUT
n = 100; % number of total items
k = 5; % number of defective items (from prevalence estimation)

%% Allocation
L_e = 3:12; % edge length of Hypercube (recomended: 3-12)
distribute_operation = "random"; % dummy yitems distribution linear, equal or random

S = 1000; % number of itterations for Monte Carlo simulation
Mean_test_reduction = [];
M_all = struct;
for L = L_e
    Test_reduction = nan(S,1);
    false_pos = nan(S,1);
    P_sus_all = nan(S,n);
    Inf_idx_all = nan(S,k);
    Num_neg_groups = nan(S,1);
    Num_pos_groups = nan(S,1);
    
    D = ceil(log(n)/log(L)); % dimension of hypercube
    %% Design Matrix
    for s = 1:S
        infected_idx = sort(randperm(n,k)); 
        Inf_idx_all(s,:) = infected_idx;
        C_label = nan(repmat(L,1,D));
        switch distribute_operation
            case "linear"
                C_label(1:n) = 1:n; % 1:N linear index
                if L^(D-1) == n
                    D = D-1;
                    C_label = nan(repmat(L,1,D));
                    C_label(1:n) = 1:n; % 1:N linear index
                end
            case "equal"
                C_label(round(linspace(1,L^D,n))) = 1:n; % equally spaced idx
                if L^(D-1) == n
                    D = D-1;
                    C_label = nan(repmat(L,1,D));
                    C_label(round(linspace(1,L^D,n))) = 1:n; % equally spaced idx 
                end
            case "random"
                rand_idx = sort(randperm(L^D,n));
                C_label(rand_idx) = 1:n; % randommly spaced idx
                if L^(D-1) == n
                    D = D-1;
                    C_label = nan(repmat(L,1,D));
                    rand_idx = sort(randperm(L^D,n));
                    C_label(rand_idx) = 1:n; % randommly spaced idx
                end
        end
    C_res = ismember(C_label,infected_idx);
    
    % decoder
    [Groups,Groups_pos,Groups_neg] = GroupsPosNeg(C_label,C_res,L);
    P_sus = HC_find_suspect(Groups,Groups_pos,Groups_neg,C_label);
    P_sus_all(s,1:numel(P_sus)) = P_sus;
    Test_reduction(s) = (1-((numel(P_sus)+(D*L))/n))*100;
    Num_neg_groups(s) = numel(Groups_neg);
    Num_pos_groups(s) = numel(Groups_pos);
    false_pos(s) = numel(P_sus)-k;
    end
    
    Pool_size = L^(D-1);
    Pool_num = D*L;
    Prevalence = (k/n)*100;
    
    Mean_test_reduction = [Mean_test_reduction, mean(Test_reduction)];
    M = zeros(L*D,n);
    for g = 1:(L*D)
        M(g,rmmissing(Groups(g,:))) = 1; % design matrix
    end   
    M_all.(strcat("L_",num2str(L))) = M;
end
%% Results
% expected test-reduction
figure("Position",[100 100 1100 700])
subplot(2,2,1)
plot(L_e,Mean_test_reduction); hold on
[~,idx] = max(Mean_test_reduction);
opt_L = idx+min(L_e)-1; % optimal edge length
xline(opt_L)
ylabel("Mean test-reduction",'Interpreter','latex','FontSize',13)
xlabel("Edge nodes L",'Interpreter','latex','FontSize',13)
legend("Mean test-reduction","Optimal edge length","Location","northwest",'Interpreter','latex','FontSize',13)

subplot(2,2,2)
histogram(Test_reduction,numel(unique(Test_reduction)),'FaceColor',[0.7 0.7 0.7]); hold on
xline(mean(Test_reduction)+std(Test_reduction),'--','Color','k','LineWidth',2)
xline(mean(Test_reduction)-std(Test_reduction),'--','Color','k','LineWidth',2)
xlabel("Test-reduction in \%",'Interpreter','latex','FontSize',13)
ylabel("Number of simulations",'Interpreter','latex','FontSize',13)
legend("Mean test-reduction likelyhood","Standard deviation","Location","northwest",'Interpreter','latex','FontSize',13)

subplot(2,2,[3 4])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
imagesc(M_all.(strcat("L_",num2str(opt_L))))
cmap = [1 1 1 %// light gray
        0  0  0]; %// white
colormap(cmap)
colorbar('Ytick',[.25 .75],'Yticklabel',[0 1]) %// only two values in colorbar
set(gca,'DataAspectRatio',[1 1 1])
ylabel('tests','interpreter','Latex','FontSize',12)
xlabel('items','interpreter','Latex','FontSize',12)
xline(1.5:(n-0.5),'alpha',0.2)
yline(1.5:(L*D-0.5),'alpha',0.2)
xticks([1,L:L:n])
yticks([1,L:L:(L*D)])
title('test matrix','Interpreter','latex','FontSize',13)

% print design matrix ~.csv
% writetable(array2table(M_all.(strcat("L_",num2str(opt_L)))),"Design_Matrix.csv")

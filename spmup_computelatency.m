function latency = spmup_computelatency(beta1,beta2)

% implement Henson et al 2002
% latency=2*1.78/(1 + exp(3.1*beta2/beta1))- 1.78

xBF.dt = 0.5;
xBF.name = 'hrf (with time and dispersion derivatives)';
xBF.length = 20;
[xBF] = spm_get_bf(xBF);

latency=2*1.78/(1 + exp(3.1.*beta2./beta1))- 1.78;
figure; subplot(2,2,4); plot(beta1*xBF.bf(:,1),'LineWidth',2);
hold on; plot(beta2*xBF.bf(:,3),'y','LineWidth',2);
plot(beta1*xBF.bf(:,1)+beta2*xBF.bf(:,3),'--k','LineWidth',2);
title('1*hrf + -1*2nd derivative','FontSize',13');
axis tight; box on; grid on


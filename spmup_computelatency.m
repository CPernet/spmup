function latency = spmup_computelatency(beta1,beta2)

% implement Henson et al 2002
% latency=2*1.78/(1 + exp(3.1*beta2/beta1))- 1.78

xBF.dt = 0.5;
xBF.name = 'hrf (with time derivative)';
xBF.length = 20;
[xBF] = spm_get_bf(xBF)

latency=2*1.78/(1 + exp(3.1.*beta2./beta1))- 1.78;
% figure; plot(beta1*xBF.bf(:,1),'--','LineWidth',2);
% hold on; plot(beta2*xBF.bf(:,2),'g--','LineWidth',2);
% plot(beta1*xBF.bf(:,1)+beta2*xBF.bf(:,2),'r','LineWidth',2);
% title(['Beta1 +, Beta2 +, latency=' num2str(latency)],'FontSize',13');
% axis tight; box on; grid on


% Script to simulate data for the parameter and model recovery analysis
% You need to have the toolbox Palamedes (https://www.palamedestoolbox.org/) on path to run this script.
% Author: Arthur S. Courtin  
% License: MIT (see LICENSE file) 

% Create function
function j=inv_logit(i)
j=1./(1+exp(-i));
end

%%  Set simulation parameters
age=[20:80 55:80];
ra=age-min(age);
h=[zeros(1,61) ones(1,26)];
n=15;

%%  Simulate data for model 1
for dataset=1:n

    clearvars -except age ra h n dataset
    %Sample hyperparameters
    %Kappa_0
    mu(1) = normrnd(0.67,0.09*4,1);
    mu(2) = normrnd(4.33,0.31*4,1);
    mu(3) = normrnd(22.48,1.18*4,1);
    mu(4) = normrnd(15.11,0.42*4,1);
    mu(5) = normrnd(1.19,0.22*4,1);
    mu(6) = normrnd(-0.23,0.15*4,1);
    mu(7) = normrnd(-1.60,0.11*4,1);
    mu(8) = normrnd(-0.22,0.15*4,1);
    mu(9) = normrnd(-3.33,0.40,1);
    mu(10) = normrnd(-3.84,0.36,1);
    %Kappa_s
    mu(11) = normrnd(0,0.54,1);
    mu(12) = normrnd(0,1.90,1);
    mu(13) = normrnd(0,7.39,1);
    mu(14) = normrnd(0,2.67,1);
    mu(15) = normrnd(0,1.11,1);
    mu(16) = normrnd(0,0.73,1);
    mu(17) = normrnd(0,0.50,1);
    mu(18) = normrnd(0,0.79,1);
    %Kappa_a
    mu(19) = normrnd(0,0.54/60,1);
    mu(20) = normrnd(0,1.90/60,1);
    mu(21) = normrnd(0,7.39/60,1);
    mu(22) = normrnd(0,2.67/60,1);
    mu(23) = normrnd(0,1.11/60,1);
    mu(24) = normrnd(0,0.73/60,1);
    mu(25) = normrnd(0,0.50/60,1);
    mu(26) = normrnd(0,0.79/60,1);
    
    %Kappa_0
    tau(1) = normrnd(0.54,0.09*4,1);
    tau(2) = normrnd(1.90,0.24*4,1);
    tau(3) = normrnd(7.39,0.83*4,1);
    tau(4) = normrnd(2.67,0.30*4,1);
    tau(5) = normrnd(1.11,0.23*4,1);
    tau(6) = normrnd(0.73,0.13*4,1);
    tau(7) = normrnd(0.50,0.13*4,1);
    tau(8) = normrnd(0.79,0.14*4,1);
    tau(9) = normrnd(0.66,0.35,1);
    tau(10) = normrnd(1.13,0.45,1);

    tau=abs(tau);
    
    %Sample participant paramaters
    for idx=1:10
        kappa_0(idx,:) = normrnd(mu(idx),tau(idx),length(age),1);
    end

    %Compute participant PF parameters
    alpha(1,:) = (kappa_0(1,:) + mu(11) .* h + mu(19) .*ra);
    alpha(2,:) = (kappa_0(2,:) + mu(12) .* h + mu(20) .*ra);
    alpha(3,:) = (kappa_0(3,:) + mu(13) .* h + mu(21) .*ra);
    alpha(4,:) = (kappa_0(4,:) + mu(14) .* h + mu(22) .*ra);
    beta(1,:) = exp(kappa_0(5,:) + mu(15) .* h + mu(23) .*ra);
    beta(2,:) = exp(kappa_0(6,:) + mu(16) .* h + mu(24) .*ra);
    beta(3,:) = exp(kappa_0(7,:) + mu(17) .* h + mu(25) .*ra);
    beta(4,:) = exp(kappa_0(8,:) + mu(18) .* h + mu(26) .*ra);
    gamma = .5*inv_logit(kappa_0(9,:));
    lambda = .5*inv_logit(kappa_0(10,:));
         
    save(['parameter_recovery/simulated_runs/G_uc/dataset_information_',num2str(dataset),'.mat'],"mu","tau","kappa_0","alpha","beta","gamma","lambda","age","h")

    %Simulate trial level data
    for ppt=1:length(age)
        tic
        fprintf('Dataset %i - Participant %i \n',dataset,ppt)
        PM=PAL_AMPM_setupPM( ...
            'priorAlphaRange',0:0.5:30, ...
            'priorBetaRange',linspace(log10(0.001),log10(6),51), ...
            'priorGammaRange',0.05, ...
            'priorLambdaRange',0.05, ...
            'stimRange',0:0.5:30,...
            'PF',@PAL_CumulativeNormal);
        for trial=1:30
            x=PM.xCurrent;
            y=binornd(1,gamma(ppt)+(1-gamma(ppt)-lambda(ppt))*normcdf(x,alpha(1,ppt),1/beta(1,ppt)));
            PM=PAL_AMPM_updatePM(PM,y);
        end
        save(['parameter_recovery/simulated_runs/G_uc/',num2str(dataset),'_',num2str(ppt),'_CDT.mat'],"PM")

        PM=PAL_AMPM_setupPM( ...
            'priorAlphaRange',0:0.5:20, ...
            'priorBetaRange',linspace(log10(0.001),log10(6),51), ...
            'priorGammaRange',0.05, ...
            'priorLambdaRange',0.05, ...
            'stimRange',0:0.5:20,...
            'PF',@PAL_CumulativeNormal);
        for trial=1:30
            x=PM.xCurrent;
            y=binornd(1,gamma(ppt)+(1-gamma(ppt)-lambda(ppt))*normcdf(x,alpha(2,ppt),1/beta(2,ppt)));
            PM=PAL_AMPM_updatePM(PM,y);
        end
        save(['parameter_recovery/simulated_runs/G_uc/',num2str(dataset),'_',num2str(ppt),'_WDT.mat'],"PM")

        PM=PAL_AMPM_setupPM( ...
            'priorAlphaRange',0:1:30, ...
            'priorBetaRange',linspace(log10(0.001),log10(6),51), ...
            'priorGammaRange',0.05, ...
            'priorLambdaRange',0.05, ...
            'stimRange',0:1:30,...
            'PF',@PAL_CumulativeNormal);
        for trial=1:40
            x=PM.xCurrent;
            y=binornd(1,gamma(ppt)+(1-gamma(ppt)-lambda(ppt))*normcdf(x,alpha(3,ppt),1/beta(3,ppt)));
            PM=PAL_AMPM_updatePM(PM,y);
        end
        save(['parameter_recovery/simulated_runs/G_uc/',num2str(dataset),'_',num2str(ppt),'_CPT.mat'],"PM")

        PM=PAL_AMPM_setupPM( ...
            'priorAlphaRange',0:1:20, ...
            'priorBetaRange',linspace(log10(0.001),log10(6),51), ...
            'priorGammaRange',0.05, ...
            'priorLambdaRange',0.05, ...
            'stimRange',0:1:20,...
            'PF',@PAL_CumulativeNormal);
        for trial=1:40
            x=PM.xCurrent;
            y=binornd(1,gamma(ppt)+(1-gamma(ppt)-lambda(ppt))*normcdf(x,alpha(4,ppt),1/beta(4,ppt)));
            PM=PAL_AMPM_updatePM(PM,y);
        end
        save(['parameter_recovery/simulated_runs/G_uc/',num2str(dataset),'_',num2str(ppt),'_HPT.mat'],"PM")
        toc
    end
end

%%  Simulate data for model 2
for dataset=1:n

    clearvars -except age ra h n dataset
    
    %Sample hyperparameters
    %Kappa_0
    mu(1) = normrnd(-0.52,0.15*4,1);
    mu(2) = normrnd(1.35,0.06*4,1);
    mu(3) = normrnd(3.06,0.08*4,1);
    mu(4) = normrnd(2.71,0.03*4,1);
    mu(5) = normrnd(1.20,0.20*4,1);
    mu(6) = normrnd(-0.24,0.15*4,1);
    mu(7) = normrnd(-1.64,0.12*4,1);
    mu(8) = normrnd(-0.30,0.16*4,1);
    mu(9) = normrnd(-3.32,0.41,1);
    mu(10) = normrnd(-3.74,0.35,1);
    %Kappa_s
    mu(11) = normrnd(0,0.92,1);
    mu(12) = normrnd(0,0.53,1);
    mu(13) = normrnd(0,0.38,1);
    mu(14) = normrnd(0,0.17,1);
    mu(15) = normrnd(0,1.09,1);
    mu(16) = normrnd(0,0.75,1);
    mu(17) = normrnd(0,0.53,1);
    mu(18) = normrnd(0,0.85,1);
    %Kappa_a
    mu(19) = normrnd(0,0.92/60,1);
    mu(20) = normrnd(0,0.53/60,1);
    mu(21) = normrnd(0,0.38/60,1);
    mu(22) = normrnd(0,0.17/60,1);
    mu(23) = normrnd(0,1.09/60,1);
    mu(24) = normrnd(0,0.75/60,1);
    mu(25) = normrnd(0,0.53/60,1);
    mu(26) = normrnd(0,0.85/60,1);
    
    %Kappa_0
    tau(1) = normrnd(0.92,0.13*4,1);
    tau(2) = normrnd(0.53,0.07*4,1);
    tau(3) = normrnd(0.38,0.04*4,1);
    tau(4) = normrnd(0.17,0.04*4,1);
    tau(5) = normrnd(1.09,0.19*4,1);
    tau(6) = normrnd(0.75,0.14*4,1);
    tau(7) = normrnd(0.53,0.11*4,1);
    tau(8) = normrnd(0.85,0.14*4,1);
    tau(9) = normrnd(0.63,0.36,1);
    tau(10) = normrnd(0.55,0.38,1);

    tau=abs(tau);

    %Sample participant paramaters
    for idx=1:10
        kappa_0(idx,:) = normrnd(mu(idx),tau(idx),length(age),1);
    end

    %Compute participant PF parameters
    alpha(1,:) = exp(kappa_0(1,:) + mu(11) .* h + mu(19) .*ra);
    alpha(2,:) = exp(kappa_0(2,:) + mu(12) .* h + mu(20) .*ra);
    alpha(3,:) = exp(kappa_0(3,:) + mu(13) .* h + mu(21) .*ra);
    alpha(4,:) = exp(kappa_0(4,:) + mu(14) .* h + mu(22) .*ra);
    beta(1,:) = exp(kappa_0(5,:) + mu(15) .* h + mu(23) .*ra);
    beta(2,:) = exp(kappa_0(6,:) + mu(16) .* h + mu(24) .*ra);
    beta(3,:) = exp(kappa_0(7,:) + mu(17) .* h + mu(25) .*ra);
    beta(4,:) = exp(kappa_0(8,:) + mu(18) .* h + mu(26) .*ra);
    gamma = .5*inv_logit(kappa_0(9,:));
    lambda = .5*inv_logit(kappa_0(10,:));
        
    save(['parameter_recovery/simulated_runs/G_c/dataset_information_',num2str(dataset),'.mat'],"mu","tau","kappa_0","alpha","beta","gamma","lambda","age","h")
    
    %Simulate trial level data
    for ppt=1:length(age)
        tic
        fprintf('Dataset %i - Participant %i \n',dataset,ppt)
        PM=PAL_AMPM_setupPM( ...
            'priorAlphaRange',0:0.5:30, ...
            'priorBetaRange',linspace(log10(0.001),log10(6),51), ...
            'priorGammaRange',0.05, ...
            'priorLambdaRange',0.05, ...
            'stimRange',0:0.5:30,...
            'PF',@PAL_CumulativeNormal);
        for trial=1:30
            x=PM.xCurrent;
            y=binornd(1,gamma(ppt)+(1-gamma(ppt)-lambda(ppt))*normcdf(x,alpha(1,ppt),1/beta(1,ppt)));
            PM=PAL_AMPM_updatePM(PM,y);
        end
        save(['parameter_recovery/simulated_runs/G_c/',num2str(dataset),'_',num2str(ppt),'_CDT.mat'],"PM")

        PM=PAL_AMPM_setupPM( ...
            'priorAlphaRange',0:0.5:20, ...
            'priorBetaRange',linspace(log10(0.001),log10(6),51), ...
            'priorGammaRange',0.05, ...
            'priorLambdaRange',0.05, ...
            'stimRange',0:0.5:20,...
            'PF',@PAL_CumulativeNormal);
        for trial=1:30
            x=PM.xCurrent;
            y=binornd(1,gamma(ppt)+(1-gamma(ppt)-lambda(ppt))*normcdf(x,alpha(2,ppt),1/beta(2,ppt)));
            PM=PAL_AMPM_updatePM(PM,y);
        end
        save(['parameter_recovery/simulated_runs/G_c/',num2str(dataset),'_',num2str(ppt),'_WDT.mat'],"PM")

        PM=PAL_AMPM_setupPM( ...
            'priorAlphaRange',0:1:30, ...
            'priorBetaRange',linspace(log10(0.001),log10(6),51), ...
            'priorGammaRange',0.05, ...
            'priorLambdaRange',0.05, ...
            'stimRange',0:1:30,...
            'PF',@PAL_CumulativeNormal);
        for trial=1:40
            x=PM.xCurrent;
            y=binornd(1,gamma(ppt)+(1-gamma(ppt)-lambda(ppt))*normcdf(x,alpha(3,ppt),1/beta(3,ppt)));
            PM=PAL_AMPM_updatePM(PM,y);
        end
        save(['parameter_recovery/simulated_runs/G_c/',num2str(dataset),'_',num2str(ppt),'_CPT.mat'],"PM")

        PM=PAL_AMPM_setupPM( ...
            'priorAlphaRange',0:1:20, ...
            'priorBetaRange',linspace(log10(0.001),log10(6),51), ...
            'priorGammaRange',0.05, ...
            'priorLambdaRange',0.05, ...
            'stimRange',0:1:20,...
            'PF',@PAL_CumulativeNormal);
        for trial=1:40
            x=PM.xCurrent;
            y=binornd(1,gamma(ppt)+(1-gamma(ppt)-lambda(ppt))*normcdf(x,alpha(4,ppt),1/beta(4,ppt)));
            PM=PAL_AMPM_updatePM(PM,y);
        end
        save(['parameter_recovery/simulated_runs/G_c/',num2str(dataset),'_',num2str(ppt),'_HPT.mat'],"PM")
        toc
    end
end

%%  Simulate data for model 3
for dataset=1:n

    clearvars -except age ra h n dataset
    
    %Sample hyperparameters
    %Kappa_0
    mu(1) = normrnd(-0.57,0.17*4,1);
    mu(2) = normrnd(1.35,0.09*4,1);
    mu(3) = normrnd(3.09,0.06*4,1);
    mu(4) = normrnd(2.72,0.03*4,1);
    mu(5) = normrnd(0.71,0.14*4,1);
    mu(6) = normrnd(1.27,0.15*4,1);
    mu(7) = normrnd(1.51,0.11*4,1);
    mu(8) = normrnd(2.59,0.17*4,1);
    mu(9) = normrnd(-3.25,0.39,1);
    mu(10) = normrnd(-4.02,0.36,1);
    %Kappa_s
    mu(11) = normrnd(0,1.02,1);
    mu(12) = normrnd(0,0.53,1);
    mu(13) = normrnd(0,0.41,1);
    mu(14) = normrnd(0,0.16,1);
    mu(15) = normrnd(0,0.39,1);
    mu(16) = normrnd(0,0.73,1);
    mu(17) = normrnd(0,0.33,1);
    mu(18) = normrnd(0,0.97,1);
    %Kappa_a
    mu(19) = normrnd(0,1.02/60,1);
    mu(20) = normrnd(0,0.53/60,1);
    mu(21) = normrnd(0,0.41/60,1);
    mu(22) = normrnd(0,0.16/60,1);
    mu(23) = normrnd(0,0.39/60,1);
    mu(24) = normrnd(0,0.73/60,1);
    mu(25) = normrnd(0,0.33/60,1);
    mu(26) = normrnd(0,0.97/60,1);
    
    %Kappa_0
    tau(1) = normrnd(1.02,0.15*4,1);
    tau(2) = normrnd(0.53,0.07*4,1);
    tau(3) = normrnd(0.41,0.05*4,1);
    tau(4) = normrnd(0.16,0.04*4,1);
    tau(5) = normrnd(0.39,0.19*4,1);
    tau(6) = normrnd(0.73,0.14*4,1);
    tau(7) = normrnd(0.33,0.14*4,1);
    tau(8) = normrnd(0.97,0.15*4,1);
    tau(9) = normrnd(0.86,0.36,1);
    tau(10) = normrnd(0.39,0.28,1);

    tau=abs(tau);

    %Sample participant paramaters
    for idx=1:10
        kappa_0(idx,:) = normrnd(mu(idx),tau(idx),length(age),1);
    end
    
    %Compute participant PF parameters
    alpha(1,:) = exp(kappa_0(1,:) + mu(11) .* h + mu(19) .*ra);
    alpha(2,:) = exp(kappa_0(2,:) + mu(12) .* h + mu(20) .*ra);
    alpha(3,:) = exp(kappa_0(3,:) + mu(13) .* h + mu(21) .*ra);
    alpha(4,:) = exp(kappa_0(4,:) + mu(14) .* h + mu(22) .*ra);
    beta(1,:) = exp(kappa_0(5,:) + mu(15) .* h + mu(23) .*ra);
    beta(2,:) = exp(kappa_0(6,:) + mu(16) .* h + mu(24) .*ra);
    beta(3,:) = exp(kappa_0(7,:) + mu(17) .* h + mu(25) .*ra);
    beta(4,:) = exp(kappa_0(8,:) + mu(18) .* h + mu(26) .*ra);
    gamma = .5*inv_logit(kappa_0(9,:));
    lambda = .5*inv_logit(kappa_0(10,:));
        
    save(['parameter_recovery/simulated_runs/Q_c/dataset_information_',num2str(dataset),'.mat'],"mu","tau","kappa_0","alpha","beta","gamma","lambda","age","h")

    %Simulate trial level data
    for ppt=1:length(age)
        tic
        fprintf('Dataset %i - Participant %i \n',dataset,ppt)
        PM=PAL_AMPM_setupPM( ...
            'priorAlphaRange',0:0.5:30, ...
            'priorBetaRange',linspace(log10(0.001),log10(6),51), ...
            'priorGammaRange',0.05, ...
            'priorLambdaRange',0.05, ...
            'stimRange',0:0.5:30,...
            'PF',@PAL_CumulativeNormal);
        for trial=1:30
            x=PM.xCurrent;
            y=binornd(1,gamma(ppt)+(1-gamma(ppt)-lambda(ppt))*(1-2^(-(x./alpha(1,ppt))^beta(1,ppt))));
            PM=PAL_AMPM_updatePM(PM,y);
        end
        save(['parameter_recovery/simulated_runs/Q_c/',num2str(dataset),'_',num2str(ppt),'_CDT.mat'],"PM")

        PM=PAL_AMPM_setupPM( ...
            'priorAlphaRange',0:0.5:20, ...
            'priorBetaRange',linspace(log10(0.001),log10(6),51), ...
            'priorGammaRange',0.05, ...
            'priorLambdaRange',0.05, ...
            'stimRange',0:0.5:20,...
            'PF',@PAL_CumulativeNormal);
        for trial=1:30
            x=PM.xCurrent;
            y=binornd(1,gamma(ppt)+(1-gamma(ppt)-lambda(ppt))*(1-2^(-(x./alpha(2,ppt))^beta(2,ppt))));
            PM=PAL_AMPM_updatePM(PM,y);
        end
        save(['parameter_recovery/simulated_runs/Q_c/',num2str(dataset),'_',num2str(ppt),'_WDT.mat'],"PM")

        PM=PAL_AMPM_setupPM( ...
            'priorAlphaRange',0:1:30, ...
            'priorBetaRange',linspace(log10(0.001),log10(6),51), ...
            'priorGammaRange',0.05, ...
            'priorLambdaRange',0.05, ...
            'stimRange',0:1:30,...
            'PF',@PAL_CumulativeNormal);
        for trial=1:40
            x=PM.xCurrent;
            y=binornd(1,gamma(ppt)+(1-gamma(ppt)-lambda(ppt))*(1-2^(-(x./alpha(3,ppt))^beta(3,ppt))));
            PM=PAL_AMPM_updatePM(PM,y);
        end
        save(['parameter_recovery/simulated_runs/Q_c/',num2str(dataset),'_',num2str(ppt),'_CPT.mat'],"PM")

        PM=PAL_AMPM_setupPM( ...
            'priorAlphaRange',0:1:20, ...
            'priorBetaRange',linspace(log10(0.001),log10(6),51), ...
            'priorGammaRange',0.05, ...
            'priorLambdaRange',0.05, ...
            'stimRange',0:1:20,...
            'PF',@PAL_CumulativeNormal);
        for trial=1:40
            x=PM.xCurrent;
            y=binornd(1,gamma(ppt)+(1-gamma(ppt)-lambda(ppt))*(1-2^(-(x./alpha(4,ppt))^beta(4,ppt))));
            PM=PAL_AMPM_updatePM(PM,y);
        end
        save(['parameter_recovery/simulated_runs/Q_c/',num2str(dataset),'_',num2str(ppt),'_HPT.mat'],"PM")
        toc
    end
end
% Script to aggregate the data simulated for parameter and model recovery analysis
% Author: Arthur S. Courtin  
% License: MIT (see LICENSE file) 

%% Set simulation parameters
close all
clear all

age=[20:80 55:80];
ra=age-min(age);
h=[zeros(1,61) ones(1,26)];
n=15;

%%  Model 1
Guc_dataset=[];
Guc_participant=[];
Guc_trial=[];

%loop through simualted datasets
for dataset=1:n
    %load dataset
    load(['parameter_recovery/simulated_runs/G_uc/dataset_information_',num2str(dataset),'.mat'])

    %append dataset dataset information
    Guc_dataset=[...
        Guc_dataset;...
        dataset mu tau...
        ];

    %append participant information
    Guc_participant=[...
        Guc_participant; ...
        repmat(dataset,size(kappa_0,2),1) (1:length(age))'  age' h' repmat(1,size(kappa_0,2),1) kappa_0(1,:)' alpha(1,:)' beta(1,:)' gamma(:) lambda(:);...
        repmat(dataset,size(kappa_0,2),1) (1:length(age))'  age' h' repmat(2,size(kappa_0,2),1) kappa_0(2,:)' alpha(2,:)' beta(2,:)' gamma(:) lambda(:);...
        repmat(dataset,size(kappa_0,2),1) (1:length(age))'  age' h' repmat(3,size(kappa_0,2),1) kappa_0(3,:)' alpha(3,:)' beta(3,:)' gamma(:) lambda(:);...
        repmat(dataset,size(kappa_0,2),1) (1:length(age))'  age' h' repmat(4,size(kappa_0,2),1) kappa_0(4,:)' alpha(4,:)' beta(4,:)' gamma(:) lambda(:);...
        ];

    %load and append trial information
    for ppt=1:length(age)
        tic
        fprintf('Dataset %i - Participant %i \n',dataset,ppt)
        
        load(['parameter_recovery/simulated_runs/G_uc/',num2str(dataset),'_',num2str(ppt),'_CDT.mat'])
        Guc_trial=[...
            Guc_trial;...
            repmat(dataset,30,1) repmat(ppt,30,1) repmat(1,30,1) PM.x(1:end-1)' PM.response'];

        load(['parameter_recovery/simulated_runs/G_uc/',num2str(dataset),'_',num2str(ppt),'_WDT.mat'])
        Guc_trial=[...
            Guc_trial;...
            repmat(dataset,30,1) repmat(ppt,30,1) repmat(2,30,1) PM.x(1:end-1)' PM.response'];
         
        load(['parameter_recovery/simulated_runs/G_uc/',num2str(dataset),'_',num2str(ppt),'_CPT.mat'])
        Guc_trial=[...
            Guc_trial;...
            repmat(dataset,40,1) repmat(ppt,40,1) repmat(3,40,1) PM.x(1:end-1)' PM.response'];
        
        load(['parameter_recovery/simulated_runs/G_uc/',num2str(dataset),'_',num2str(ppt),'_HPT.mat'])
        Guc_trial=[...
            Guc_trial;...
            repmat(dataset,40,1) repmat(ppt,40,1) repmat(4,40,1) PM.x(1:end-1)' PM.response'];

        toc
    end
end

%Save aggregated data
Guc_dataset=array2table(Guc_dataset);
writetable(Guc_dataset,'parameter_recovery/simulated_runs/Guc_dataset_info.csv');

Guc_participant=array2table(Guc_participant,'VariableNames',{'dataset','participant','age','status','type','kappa_0','alpha','beta','gamma','lambda'});
writetable(Guc_participant,'parameter_recovery/simulated_runs/Guc_participant_info.csv');

Guc_trial=array2table(Guc_trial,'VariableNames',{'dataset','participant','type','x','response'});
writetable(Guc_trial,'parameter_recovery/simulated_runs/Guc_trial_info.csv');

%%  Model 2
Gc_dataset=[];
Gc_participant=[];
Gc_trial=[];

%loop through simualted datasets
for dataset=1:n
    %load dataset
    load(['parameter_recovery/simulated_runs/G_c/dataset_information_',num2str(dataset),'.mat'])

    %append dataset dataset information
    Gc_dataset=[...
        Gc_dataset;...
        dataset mu tau...
        ];

    %append participant information
    Gc_participant=[...
        Gc_participant; ...
        repmat(dataset,size(kappa_0,2),1) (1:length(age))'  age' h' repmat(1,size(kappa_0,2),1) kappa_0(1,:)' alpha(1,:)' beta(1,:)' gamma(:) lambda(:);...
        repmat(dataset,size(kappa_0,2),1) (1:length(age))'  age' h' repmat(2,size(kappa_0,2),1) kappa_0(2,:)' alpha(2,:)' beta(2,:)' gamma(:) lambda(:);...
        repmat(dataset,size(kappa_0,2),1) (1:length(age))'  age' h' repmat(3,size(kappa_0,2),1) kappa_0(3,:)' alpha(3,:)' beta(3,:)' gamma(:) lambda(:);...
        repmat(dataset,size(kappa_0,2),1) (1:length(age))'  age' h' repmat(4,size(kappa_0,2),1) kappa_0(4,:)' alpha(4,:)' beta(4,:)' gamma(:) lambda(:);...
        ];

    %load and append trial information
    for ppt=1:length(age)
        tic
        fprintf('Dataset %i - Participant %i \n',dataset,ppt)
        
        load(['parameter_recovery/simulated_runs/G_c/',num2str(dataset),'_',num2str(ppt),'_CDT.mat'])
        Gc_trial=[...
            Gc_trial;...
            repmat(dataset,30,1) repmat(ppt,30,1) repmat(1,30,1) PM.x(1:end-1)' PM.response'];

        load(['parameter_recovery/simulated_runs/G_c/',num2str(dataset),'_',num2str(ppt),'_WDT.mat'])
        Gc_trial=[...
            Gc_trial;...
            repmat(dataset,30,1) repmat(ppt,30,1) repmat(2,30,1) PM.x(1:end-1)' PM.response'];
         
        load(['parameter_recovery/simulated_runs/G_c/',num2str(dataset),'_',num2str(ppt),'_CPT.mat'])
        Gc_trial=[...
            Gc_trial;...
            repmat(dataset,40,1) repmat(ppt,40,1) repmat(3,40,1) PM.x(1:end-1)' PM.response'];
        
        load(['parameter_recovery/simulated_runs/G_c/',num2str(dataset),'_',num2str(ppt),'_HPT.mat'])
        Gc_trial=[...
            Gc_trial;...
            repmat(dataset,40,1) repmat(ppt,40,1) repmat(4,40,1) PM.x(1:end-1)' PM.response'];

        toc
    end
end

%Save aggregated data
Gc_dataset=array2table(Gc_dataset);
writetable(Gc_dataset,'parameter_recovery/simulated_runs/Gc_dataset_info.csv');

Gc_participant=array2table(Gc_participant,'VariableNames',{'dataset','participant','age','status','type','kappa_0','alpha','beta','gamma','lambda'});
writetable(Gc_participant,'parameter_recovery/simulated_runs/Gc_participant_info.csv');

Gc_trial=array2table(Gc_trial,'VariableNames',{'dataset','participant','type','x','response'});
writetable(Gc_trial,'parameter_recovery/simulated_runs/.csv');

%%  Model 3
Qc_dataset=[];
Qc_participant=[];
Qc_trial=[];

%loop through simualted datasets
for dataset=1:n
    %load dataset
    load(['parameter_recovery/simulated_runs/Q_c/dataset_information_',num2str(dataset),'.mat'])

    %append dataset dataset information
    Qc_dataset=[...
        Qc_dataset;...
        dataset mu tau...
        ];

    %append participant information
    Qc_participant=[...
        Qc_participant; ...
        repmat(dataset,size(kappa_0,2),1) (1:length(age))'  age' h' repmat(1,size(kappa_0,2),1) kappa_0(1,:)' alpha(1,:)' beta(1,:)' gamma(:) lambda(:);...
        repmat(dataset,size(kappa_0,2),1) (1:length(age))'  age' h' repmat(2,size(kappa_0,2),1) kappa_0(2,:)' alpha(2,:)' beta(2,:)' gamma(:) lambda(:);...
        repmat(dataset,size(kappa_0,2),1) (1:length(age))'  age' h' repmat(3,size(kappa_0,2),1) kappa_0(3,:)' alpha(3,:)' beta(3,:)' gamma(:) lambda(:);...
        repmat(dataset,size(kappa_0,2),1) (1:length(age))'  age' h' repmat(4,size(kappa_0,2),1) kappa_0(4,:)' alpha(4,:)' beta(4,:)' gamma(:) lambda(:);...
        ];

    %load and append trial information
    for ppt=1:length(age)
        tic
        fprintf('Dataset %i - Participant %i \n',dataset,ppt)
        
        load(['parameter_recovery/simulated_runs/Q_c/',num2str(dataset),'_',num2str(ppt),'_CDT.mat'])
        Qc_trial=[...
            Qc_trial;...
            repmat(dataset,30,1) repmat(ppt,30,1) repmat(1,30,1) PM.x(1:end-1)' PM.response'];

        load(['parameter_recovery/simulated_runs/Q_c/',num2str(dataset),'_',num2str(ppt),'_WDT.mat'])
        Qc_trial=[...
            Qc_trial;...
            repmat(dataset,30,1) repmat(ppt,30,1) repmat(2,30,1) PM.x(1:end-1)' PM.response'];
         
        load(['parameter_recovery/simulated_runs/Q_c/',num2str(dataset),'_',num2str(ppt),'_CPT.mat'])
        Qc_trial=[...
            Qc_trial;...
            repmat(dataset,40,1) repmat(ppt,40,1) repmat(3,40,1) PM.x(1:end-1)' PM.response'];
        
        load(['parameter_recovery/simulated_runs/Q_c/',num2str(dataset),'_',num2str(ppt),'_HPT.mat'])
        Qc_trial=[...
            Qc_trial;...
            repmat(dataset,40,1) repmat(ppt,40,1) repmat(4,40,1) PM.x(1:end-1)' PM.response'];

        toc
    end
end

%Save aggregated data
Qc_dataset=array2table(Qc_dataset);
writetable(Qc_dataset,'parameter_recovery/simulated_runs/Qc_dataset_info.csv');

Qc_participant=array2table(Qc_participant,'VariableNames',{'dataset','participant','age','status','type','kappa_0','alpha','beta','gamma','lambda'});
writetable(Qc_participant,'parameter_recovery/simulated_runs/Qc_participant_info.csv');

Qc_trial=array2table(Qc_trial,'VariableNames',{'dataset','participant','type','x','response'});
writetable(Qc_trial,'parameter_recovery/simulated_runs/Qc_trial_info.csv');

clear all
list=ls('fMRI experiment data\');
curr_dir=pwd();
for s=1:length(list)-2
    load([curr_dir '\fMRI experiment data\' list(s+2,:) '\ses-01\beh\01_thr_detect\' list(s+2,:) '_ses-01_task-psidetcold_beh.mat']);
    cd(s).x=Results.PM.x(1:end-1);
    cd(s).y=Results.PM.response;
    cd(s).temp=Results.tcsData;
    clearvars -except s cd wd cp hp curr_dir list
    
    load([curr_dir '\fMRI experiment data\' list(s+2,:) '\ses-01\beh\01_thr_detect\' list(s+2,:) '_ses-01_task-psidetwarm_beh.mat']);
    wd(s).x=Results.PM.x(1:end-1);
    wd(s).y=Results.PM.response;
    wd(s).temp=Results.tcsData;
    clearvars -except s cd wd cp hp curr_dir list

    load([curr_dir '\fMRI experiment data\' list(s+2,:) '\ses-01\beh\02_thr_pain\' list(s+2,:) '_ses-01_task-psipaincold_beh.mat']);
    cp(s).x=Results.PM.x(1:end-1);
    cp(s).y=Results.PM.response;
    cp(s).temp=Results.tcsData;
    clearvars -except s cd wd cp hp curr_dir list

    load([curr_dir '\fMRI experiment data\' list(s+2,:) '\ses-01\beh\02_thr_pain\' list(s+2,:) '_ses-01_task-psipainwarm_beh.mat']);
    hp(s).x=Results.PM.x(1:end-1);
    hp(s).y=Results.PM.response;
    hp(s).temp=Results.tcsData;
    clearvars -except s cd wd cp hp curr_dir list
end
save('aggregated_data.mat','hp','cp','wd','cd')

%%
clc
load('aggregated_data.mat')

x=[];
y=[];
p=[];
c=[];

cdm=[];
cdv=[];
wdm=[];
wdv=[];
cpm=[];
cpv=[];
hpm=[];
hpv=[];

for s=1:length(cd)
    for t=1:min(length(cd(s).x),length(cd(s).temp))
        temps=cell2mat(cd(s).temp(t));
        temps=temps(:,2:end);
        temps(:,2:end)=abs(temps(:,2:end)-30);
        temps(:,1)=temps(:,1)-temps(1,1);

        d=max(temps(:,1))-.5;
        idx=find(temps(:,1)>d,1);
        mt=mean(temps(idx:end,2:end));
        vt=mean(abs(temps(idx:end,2:end)-mt));
        cdm=[cdm maxk(mt,2)];
        cdv=[cdv mean([vt(mt==cdm(end-1)) vt(mt==cdm(end))])];
        
        if ~isnan(cd(s).y(t)) && abs(cdm(end)-abs(cd(s).x(t)-30))<1
            x=[x mean(cdm(end-1:end))];
            y=[y abs(cd(s).y(t)-1)];
            p=[p s];
            c=[c 1];
        end
    end
    for t=1:min(length(wd(s).x),length(wd(s).temp))
        temps=cell2mat(wd(s).temp(t));
        temps=temps(:,2:end);
        temps(:,2:end)=abs(temps(:,2:end)-30);
        temps(:,1)=temps(:,1)-temps(1,1);
      
        d=max(temps(:,1))-.5;
        idx=find(temps(:,1)>d,1);
        mt=mean(temps(idx:end,2:end));
        vt=mean(abs(temps(idx:end,2:end)-mt));
        wdm=[wdm maxk(mt,2)];
        wdv=[wdv mean([vt(mt==wdm(end-1)) vt(mt==wdm(end))])];
        
        if ~isnan(wd(s).y(t)) && abs(wdm(end)-abs(wd(s).x(t)-30))<1
            x=[x mean(wdm(end-1:end))];
            y=[y wd(s).y(t)];
            p=[p s];
            c=[c 2];
        end
    end
    for t=1:min(length(cp(s).x),length(cp(s).temp))
        temps=cell2mat(cp(s).temp(t+10));
        temps=temps(:,2:end);
        temps(:,2:end)=abs(temps(:,2:end)-30);
        temps(:,1)=temps(:,1)-temps(1,1);

        d=max(temps(:,1))-.5;
        idx=find(temps(:,1)>d,1);
        mt=mean(temps(idx:end,2:end));
        vt=mean(abs(temps(idx:end,2:end)-mt));
        cpm=[cpm maxk(mt,2)];
        cpv=[cpv mean([vt(mt==cpm(end-1)) vt(mt==cpm(end))])];
        
        if ~isnan(cp(s).y(t)) && abs(cpm(end)-abs(cp(s).x(t)-30))<1
            x=[x mean(cpm(end-1:end))];
            y=[y abs(cp(s).y(t)-1)];
            p=[p s];
            c=[c 3];
        end
    end
    for t=1:min(length(hp(s).x),length(hp(s).temp))
        temps=cell2mat(hp(s).temp(t+10));
        temps=temps(:,2:end);
        temps(:,2:end)=abs(temps(:,2:end)-30);
        temps(:,1)=temps(:,1)-temps(1,1);
       
        d=max(temps(:,1))-.5;
        idx=find(temps(:,1)>d,1);
        mt=mean(temps(idx:end,2:end));
        vt=mean(abs(temps(idx:end,2:end)-mt));
        hpm=[hpm maxk(mt,2)];
        hpv=[hpv mean([vt(mt==hpm(end-1)) vt(mt==hpm(end))])];
        
        if ~isnan(hp(s).y(t)) && abs(hpm(end)-abs(hp(s).x(t)-30))<1
            x=[x mean(hpm(end-1:end))];
            y=[y hp(s).y(t)];
            p=[p s];
            c=[c 4];
        end
    end
end

results_table=table(x',y',p',c','VariableNames',{'x','y','p','c'});
writetable(results_table,'alex_fmri_thresholding_data.csv')
%%
subplot(2,2,1)
histogram(cdv)
subplot(2,2,2)
histogram(wdv)
subplot(2,2,3)
histogram(cpv)
subplot(2,2,4)
histogram(hpv)

quantile(cdv,.95)
quantile(wdv,.95)
quantile(cpv,.95)
quantile(hpv,.95)
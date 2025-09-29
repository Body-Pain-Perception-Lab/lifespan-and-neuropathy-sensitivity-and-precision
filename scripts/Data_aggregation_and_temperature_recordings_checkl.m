% Script to aggregate the data from the experiment folders on cortex
% Author: Arthur S. Courtin  
% License: MIT (see LICENSE file) 

clc
clear all

path_to_dataset='/mnt/cortex/bpp/MINDLAB2024_BEH-ageing-neuropathy';
participant_to_exclude = readtable("/home/arthur/R/thermal_ageing_thresholding/data/participant_to_exclude.csv");


folder_list=dir(path_to_dataset);

T=[];
Y=[];
R=[];
M=[];
F1=[];
F2=[];
V1=[];
V2=[];
C=[];
S=[];
R={};

for idx=1:length(folder_list)
    if contains(folder_list(idx).name,'sub-')
        subject=folder_list(idx).name;
        disp(subject)
        subject_number=str2num(subject(5:8));
        
        if sum(subject_number==participant_to_exclude.ID)==0
            %CDT    
            t=[];
            y=[];
            r=[];
            m=[];
            f1=[];
            f2=[];
            v1=[];
            v2=[];        
    
            load([path_to_dataset,'/',subject,'/ses-01/beh/01_thr_detect/',subject,'_ses-01_task-psidetcold_beh.mat']);
    
            t=Results.targetTcold;
            y=Results.Response(:,1);
            r=Results.tcsData;
    
            for jdx=1:length(r)
                temp_rec=cell2mat(r(jdx));
                time=temp_rec(:,2);
                tdx=find(abs(time-(max(time)-1))==min(abs(time-(max(time)-1))));
                delta=temp_rec(tdx:end,3:end)-t(jdx);
                delta_m=abs(mean(delta));
                mins=mink(delta_m,2);
                if mins(1)==mins(2)
                    z_idx=[find(delta_m==mins(1))];
                else
                    z_idx=[find(delta_m==mins(1)) find(delta_m==mins(2))];
                end
                m(jdx)=mean(temp_rec(tdx:end,z_idx+2),[1:2]);
                f1(jdx)=sum(abs(temp_rec(tdx:end,z_idx+2)-m(jdx))>.2,[1 2])>0;
                f2(jdx)=sum(abs(delta(:,z_idx))>1,[1 2])>0;
                v1(jdx)=var(delta(:,z_idx(1)));
                v2(jdx)=var(delta(:,z_idx(2)));
            end
    
            C=[C; repmat(1,length(t),1)];
            S=[S; repmat(str2num(strrep(subject,'sub-','')),length(t),1)];            
            T=[T; 30-t];
            M=[M; 30-m'];
            Y=[Y; y]; 
            F1=[F1; f1'];
            F2=[F2; f2'];
            V1=[V1; v1'];
            V2=[V2; v2'];
            R=[R;r];
            
            %WDT   
            t=[];
            y=[];
            r=[];
            m=[];
            f1=[];
            f2=[];
            v1=[];
            v2=[];        
    
            load([path_to_dataset,'/',subject,'/ses-01/beh/01_thr_detect/',subject,'_ses-01_task-psidetwarm_beh.mat']);
    
            t=Results.targetTwarm;
            y=Results.Response(:,1);
            r=Results.tcsData;
    
            for jdx=1:length(r)
                temp_rec=cell2mat(r(jdx));
                time=temp_rec(:,2);
                tdx=find(abs(time-(max(time)-1))==min(abs(time-(max(time)-1))));
                delta=temp_rec(tdx:end,3:end)-t(jdx);
                delta_m=abs(mean(delta));
                mins=mink(delta_m,2);
                if mins(1)==mins(2)
                    z_idx=[find(delta_m==mins(1))];
                else
                    z_idx=[find(delta_m==mins(1)) find(delta_m==mins(2))];
                end
                m(jdx)=mean(temp_rec(tdx:end,z_idx+2),[1:2]);
                f1(jdx)=sum(abs(temp_rec(tdx:end,z_idx+2)-m(jdx))>.2,[1 2])>0;
                f2(jdx)=sum(abs(delta(:,z_idx))>1,[1 2])>0;
                v1(jdx)=var(delta(:,z_idx(1)));
                v2(jdx)=var(delta(:,z_idx(2)));
            end
    
            C=[C; repmat(2,length(t),1)];
            S=[S; repmat(str2num(strrep(subject,'sub-','')),length(t),1)];            
            T=[T; t-30];
            M=[M; m'-30];
            Y=[Y; y]; 
            F1=[F1; f1'];
            F2=[F2; f2'];
            V1=[V1; v1'];
            V2=[V2; v2'];
            R=[R;r];
            
            %CPT   
            t=[];
            y=[];
            r=[];
            m=[];
            f1=[];
            f2=[];
            v1=[];
            v2=[];        
    
            load([path_to_dataset,'/',subject,'/ses-01/beh/02_thr_pain/',subject,'_ses-01_task-psipaincold_beh.mat']);
    
            t=Results.targetTcold(11:end);
            y=Results.Response(11:end,1);
            r=Results.tcsData(11:end);
    
            for jdx=1:length(r)
                temp_rec=cell2mat(r(jdx));
                time=temp_rec(:,2);
                tdx=find(abs(time-(max(time)-1))==min(abs(time-(max(time)-1))));
                delta=temp_rec(tdx:end,3:end)-t(jdx);
                delta_m=abs(mean(delta));
                mins=mink(delta_m,2);
                if mins(1)==mins(2)
                    z_idx=[find(delta_m==mins(1))];
                else
                    z_idx=[find(delta_m==mins(1)) find(delta_m==mins(2))];
                end
                m(jdx)=mean(temp_rec(tdx:end,z_idx+2),[1:2]);
                f1(jdx)=sum(abs(temp_rec(tdx:end,z_idx+2)-m(jdx))>.2,[1 2])>0;
                f2(jdx)=sum(abs(delta(:,z_idx))>1,[1 2])>0;
                v1(jdx)=var(delta(:,z_idx(1)));
                v2(jdx)=var(delta(:,z_idx(2)));
            end
    
            C=[C; repmat(3,length(t),1)];
            S=[S; repmat(str2num(strrep(subject,'sub-','')),length(t),1)];            
            T=[T; 30-t];
            M=[M; 30-m'];
            Y=[Y; y]; 
            F1=[F1; f1'];
            F2=[F2; f2'];
            V1=[V1; v1'];
            V2=[V2; v2'];
            R=[R;r];
            
            %HPT   
            t=[];
            y=[];
            r=[];
            m=[];
            f1=[];
            f2=[];
            v1=[];
            v2=[];        
    
            load([path_to_dataset,'/',subject,'/ses-01/beh/02_thr_pain/',subject,'_ses-01_task-psipainwarm_beh.mat']);
    
            t=Results.targetTwarm(11:end);
            y=Results.Response(11:end,1);
            r=Results.tcsData(11:end);
    
            for jdx=1:length(r)
                temp_rec=cell2mat(r(jdx));
                time=temp_rec(:,2);
                tdx=find(abs(time-(max(time)-1))==min(abs(time-(max(time)-1))));
                delta=temp_rec(tdx:end,3:end)-t(jdx);
                delta_m=abs(mean(delta));
                mins=mink(delta_m,2);
                if mins(1)==mins(2)
                    z_idx=[find(delta_m==mins(1))];
                else
                    z_idx=[find(delta_m==mins(1)) find(delta_m==mins(2))];
                end
                m(jdx)=mean(temp_rec(tdx:end,z_idx+2),[1:2]);
                f1(jdx)=sum(abs(temp_rec(tdx:end,z_idx+2)-m(jdx))>.2,[1 2])>0;
                f2(jdx)=sum(abs(delta(:,z_idx))>1,[1 2])>0;
                v1(jdx)=var(delta(:,z_idx(1)));
                v2(jdx)=var(delta(:,z_idx(2)));
            end
    
            C=[C; repmat(4,length(t),1)];
            S=[S; repmat(str2num(strrep(subject,'sub-','')),length(t),1)];            
            T=[T; t-30];
            M=[M; m'-30];
            Y=[Y; y]; 
            F1=[F1; f1'];
            F2=[F2; f2'];
            V1=[V1; v1'];
            V2=[V2; v2'];
            R=[R;r];
        end
    end
end
%%
aggregated_results=table(C,S,T,M,Y,F1,F2,V1,V2,'VariableNames',{'condition','subject','target_intensity','recorded_intensity','response','recording_deviates_from_mean','recording_deviates_from_target','recording_variance_z1','recording_variance_z2'});
writetable(aggregated_results,'/home/arthur/R/thermal_ageing_thresholding/data/aggregated_results.csv')
%%
ddx=find(F1>0);
for idx=1:length(ddx)
   subplot(1,length(ddx),idx)
   temp_rec=cell2mat(R(ddx(idx)));
   time=temp_rec(:,2)-min(temp_rec(:,2));
   plot(time,temp_rec(:,3:end))
   title(num2str(ddx(idx)))
end
%%
V=mean([V1,V2],2);
vdx=find(V>(quantile(V,0.75)+1.5*iqr(V)));

for idx=1:length(vdx)
   close all
    figure('Units','normalized','Position',[0 0 1 1])
   temp_rec=cell2mat(R(vdx(idx)));
   time=temp_rec(:,2)-min(temp_rec(:,2));
   temp=temp_rec(:,3:end);
   plot(time,temp)
   hold on
   if mean(temp-30,'all')>0
       plot([min(time) max(time)],ones(1,2)*(M(vdx(idx))+30.2),'k:')
       plot([min(time) max(time)],ones(1,2)*(M(vdx(idx))+29.8),'k:')
   else
       plot([min(time) max(time)],ones(1,2)*(30.2-M(vdx(idx))),'k:')
       plot([min(time) max(time)],ones(1,2)*(29.8-M(vdx(idx))),'k:')
   end
   hold off
   title([num2str(vdx(idx)) ' ' num2str(M(vdx(idx)))])
   pause(1)
end


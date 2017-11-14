function [RR,brady,tachy,PVC,PAC]=ECG_diagnosis(sig,fs)
%   ECG_diagnosis  distinguish R wave and find out if there are any diseases
%regarding brady,tachy,PVC and PAC
%   sig    the original ECG wave
%   fs     sampling rate
%   RR     the RR interval
%   brady  the index that looks like bradycardia 
%   tachy  the index that looks like tachycardia
%   PVC    the index that looks like Premature ventricular contraction
%   PAC    the index that looks like Premature atrium contraction
tm=(1:length(sig))/fs;
figure();
subplot(3,2,1);plot(tm(1:1800),sig(1:1800));ylabel('Raw signal')
% lowpass
B=[1 0 0 0 0 0 -2 0 0 0 0 0 1];
A=[1 -2 1];
sig_L=filter(B,A,sig);
subplot(3,2,2);plot(tm(1:1800),sig_L(1:1800));ylabel('After lowpass');
% highpass
B=[-1/32 zeros(1,15) 1 -1 zeros(1,14) 1/32];
A=[1 -1];
sig_h=filter(B,A,sig_L);
subplot(3,2,3);plot(tm(1:1800),sig_h(1:1800));ylabel('After highpass');
% difference
B=1/8*[2 1 0 -1 -2];
A=[1];
sig_d=filter(B,A,sig_h);
subplot(3,2,4);plot(tm(1:1800),sig_d(1:1800));ylabel('Derivative')
% square
sig_s=sig_d.^2;
subplot(3,2,5);plot(tm(1:1800),sig_s(1:1800));ylabel('Squaring')
% window average
N=round(0.15*fs);
B=1/N*ones(1,N);
A=[1];
sig_w=filter(B,A,sig_s);
subplot(3,2,6);plot(tm(1:1800),sig_w(1:1800));ylabel('Window Integration')
% initialization
qrs_a=[];
qrs_i=[];
qrs_candidate=[];
qrs_raw=[];
T_check=0;
THR_s=[];
% Fiducial Mark
[pks,locs] = findpeaks(sig_w,'MINPEAKDISTANCE',round(0.2*fs));
% set initial threshold for window average
SPK=max(sig_w(1:round(2*fs)));
NPK=0;
T1=0.75*NPK+0.25*SPK;
T2=0.5*T1;
% set initial threshold for filtered signal 
SPK_s=max(sig_h(1:round(2*fs)));
NPK_s=0;
T1_s=0.75*NPK_s+0.25*SPK_s;
T2_s=0.5*T1_s;
% detect QRS
for i=1:length(pks)
    % find the peak in the filtered signal
    if (locs(i)-round(0.15*fs)>=1)
        [sig_amp,sig_index]=max(sig_h(locs(i)-round(0.15*fs):locs(i)));
    else
        if (i==1)
            [sig_amp,sig_index]=max(sig_h(1:locs(i)));
        end
    end
    sig_index=sig_index+locs(i)-round(0.15*fs)-1;
    qrs_candidate=[qrs_candidate sig_index];
    % check the RR interval 
    if (length(qrs_i)>=9)
        RR_int=diff(qrs_i);
        mean_RR1=mean(RR_int(end-7:end));
        RR_int2=RR_int(RR_int>=0.92*mean_RR1);
        RR_int2=RR_int2(RR_int2<=1.16*mean_RR1);
        if (length(RR_int2)>=8) mean_RR2=mean(RR_int2(end-7:end));
        else
            mean_RR2=mean(RR_int2);
        end
    else
        mean_RR2=mean(diff(qrs_i));
    end
    if ((~isempty(qrs_i))&&((locs(i)-qrs_i(end))>=round(1.66*mean_RR2)))
        % find missed beats
        [pks_temp,locs_temp] = max(sig_w(qrs_i(end)+ round(0.200*fs):locs(i)-round(0.200*fs))); % search back and locate the max in this interval
        locs_temp = qrs_i(end)+ round(0.200*fs) + locs_temp -1; %location 
        % check if it is a signal
        if (pks_temp > T2)
           % find the peak in the filtered signal
           [amp_temp,index_temp] = max(sig_h(locs_temp-round(0.150*fs):locs_temp));
           if (amp_temp>T2_s)
               qrs_a = [qrs_a pks_temp];
               qrs_i = [qrs_i locs_temp];
               qrs_raw=[qrs_raw index_temp+locs_temp-round(0.150*fs)-1];
               % change the threshold
               SPK=0.25*pks_temp+0.75*SPK;
               SPK_s=0.25*amp_temp+0.75*SPK;
               T1=0.75*SPK+0.25*NPK;
               T2=0.5*T1;
               T1_s=0.75*SPK_s+0.25*NPK_s;
               T2_s=0.5*T2;
           end
        end
    end
    % distinguish noise and QRS
    if (pks(i)>T1) 
        if ((~isempty(qrs_i))&&(locs(i)-qrs_i(end)<=round(0.36*fs)))
            % check if it is a T wave
            Slope_qrs=max(abs(diff(sig_h(qrs_raw(end)-round(0.075*fs):qrs_raw(end)))));
            Slope=max(abs(diff(sig_h(sig_index-round(0.075*fs):sig_index))));
            if (Slope<=0.5*Slope_qrs)
                T_check=1;
            end
        end
        if (T_check==0)
            % not a T wave
            if (sig_amp>T1_s) 
                %¡¡QRS is found
                qrs_i=[qrs_i locs(i)];
                qrs_a=[qrs_a pks(i)];
                qrs_raw=[qrs_raw sig_index];
                SPK=0.125*pks(i)+0.875*SPK;
                SPK_s=0.125*sig_amp+0.875*SPK_s;
            else
                NPK=0.125*pks(i)+0.875*NPK;
                NPK_s=0.125*sig_amp+0.875*NPK_s;
            end
        else
            % This is a T wave
            NPK=0.125*pks(i)+0.875*NPK;
            NPK_s=0.125*sig_amp+0.875*NPK_s;   
        end
    else
        % This is noise
        NPK=0.125*pks(i)+0.875*NPK;
        NPK_s=0.125*sig_amp+0.875*NPK_s;
    end
    % update the threshold
    T1=0.75*NPK+0.25*SPK;
    T2=0.5*T1;
    T1_s=0.75*NPK_s+0.25*SPK_s;
    T2_s=0.5*T1_s;
    THR_s=[THR_s T1_s];
    % reset T_check
    T_check=0;
end
% show the results
figure();
plot((1:length(sig_h))/fs,sig_h);
hold on;scatter(qrs_raw/fs,sig_h(qrs_raw));
hold on;plot(qrs_candidate/fs,THR_s);
% locate Q & S
q_loc=[];
s_loc=[];
% ignore the first and last QRS complex
for i=2:length(qrs_raw)-1
        [~,q_loc_temp]=min(sig_h(qrs_raw(i)-round(0.1*fs):qrs_raw(i)));
        q_loc=[q_loc,q_loc_temp+qrs_raw(i)-round(0.1*fs)-1];
        [~,s_loc_temp]=min(sig_h(qrs_raw(i):qrs_raw(i)+round(0.075*fs)));
        s_loc=[s_loc,s_loc_temp+qrs_raw(i)-1];
end
hold on;scatter(q_loc/fs,sig_h(q_loc));
hold on;scatter(s_loc/fs,sig_h(s_loc));
qrs_width=s_loc-q_loc;
% judge if there is a bradycardia
RR_int=diff(qrs_raw);
brady=qrs_raw(RR_int>round(1.5*fs));
AR_int=filter(1/8*ones(1,8),1,RR_int);
brady=[brady,qrs_raw(AR_int>round(1.2*fs))];
% judge if there is a tachycardia
AR_int(1:7)=8*AR_int(1:7)./(1:7);
tachy=qrs_raw(AR_int<round(0.6*fs));
base=qrs_raw(9);
tachy(tachy<base)=[];
% judge if there is PVC(Premature ventricular contraction) or 
% PAC(Premature atrium contraction)
PVC=[];
PAC=[];
for i=9:length(RR_int)-1
    AR_prev=AR_int(i-1);
    RR_now=RR_int(i);
    RR_next=RR_int(i+1);
    width=qrs_width(i);
    % note that the ith qrs wave in the qrs_width vector is actually the (i+1)th in the raw signal
    if ((RR_now<0.875*AR_prev)&&(width>round(0.12*fs))&&(abs(RR_next+RR_now-2*AR_prev)<=2))
        PVC=[PVC,qrs_raw(i+1)];
    end
    if ((RR_now<0.875*AR_prev)&&(width<round(0.12*fs))&&(RR_next+RR_now-2*AR_prev<=-2))
        PAC=[PAC,qrs_raw(i+1)];
    end
end
RR=mean(RR_int);
end
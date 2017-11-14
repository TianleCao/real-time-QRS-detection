function ECG_demo(fs)
%  ECG_demo   load data, plot and diagnose it in real-time.
%  fs         sampling rate
sig=load('./data/7.txt');% change to load other files 
len = length(sig);
plt_inter = 0.2; % plot interval
n_tasks = len/(fs*plt_inter); %total tasks for the timer
t=[0];
lh=plot(t,sin(t),'-');   
t=timer(...
    'Name','MyTimer',...
    'TimerFcn',@ExecuteTask1,... 
    'ErrorFcn',@ExecuteError,...
    'Period',plt_inter,'TasksToExecute',n_tasks,...
    'ExecutionMode','fixedrate'); 
ud=struct('linehandle',lh,'count',1,'signal',sig,'fs',fs,'len',len); 
set(t,'UserData',ud);
start(t);
end
function ExecuteTask1(obj,eventdata)
     ud=obj.UserData;
     fs=ud.fs;
     len=ud.len;
     l=ud.linehandle;
     c=ud.count;
     plt_inter=obj.Period;
     plt_count=plt_inter*fs;
     t=get(l,'XData');
     y=get(l,'YData');
     ECG=ud.signal;
     t=[t,(c:c+plt_count-1)/fs];
     y=[y,ECG(c:c+plt_count-1)'];
     set(ud.linehandle,'XData',t,'YData',y); 
     drawnow; 
     if ((c/fs)<=5) axis([0 5 -3 3]);
     else
         axis([c/fs-5 c/fs -3 3]);
     end
     ud.count=ud.count+plt_count; 
     if (mod(ud.count-1,2*fs)==0) 
         rate=rate_cal(ECG(1:ud.count-1),fs);
         fprintf('%ds real-time heart rate£º%.1f/min\n',(ud.count-1)/fs,rate);
     end
     set(obj,'UserData',ud);
     if (ud.count>=len) 
         [RR,brady,tachy,PVC,PAC]=ECG_diagnosis(ECG,fs);
         fprintf('Average heart rate£º%.1f/min\n',60/(RR/fs));
         if (~isempty(brady)) 
             fprintf('There seems to be Bradycardia, %d records in all have been found \n at time=',length(brady));
             fprintf('%.1fs ',brady/fs);
         end
         if (~isempty(tachy)) 
             fprintf('There seems to be Tachycardia, %d records in all have been found \n at time=',length(tachy));
             fprintf('%.1fs ',tachy/fs);
         end
         if (~isempty(PVC)) 
             fprintf('There seems to be Premature ventricular contractions, %d records in all have been found \n at time=',length(PVC));
             fprintf('%.1fs ',PVC/fs);
         end
         if (~isempty(PAC)) 
             fprintf('There seems to be Atrial premature beats, %d records in all have been found \n at time=',length(PAC));
             fprintf('%.1fs ',PAC/fs);
         end
     end
end
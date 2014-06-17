function evalCRF(leaveOut,res)
%k = 16; %silly way to avoid testing on training
%while k == leaveOut
%    k = randi(18,1); %hard coded for now
%end
k = 15;
%maunally inputed these from the ouput of tuning
starts = [400 400 400 400 400 400 400 400 500 500 500 500 200 200 200 200 200 200] 
CRFGM_paramEval(leaveOut+1,1,starts(leaveOut + 1),leaveOut+1,res);
%fprintf('tuned to %d',iter);
%CRFGM_evalIter(leaveOut,tuneResults(leaveOut+1),res);
end

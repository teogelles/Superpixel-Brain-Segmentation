function tuneCRF(leaveOut,res)
k = leaveOut; %silly way to avoid testing on training
while k == leaveOut
    k = randi(18,1); %hard coded for now to size of data set (18 images)
end
iter = CRFGM_fastTune(k,1,1000,leaveOut,res);
fprintf('tuned to %d',iter);
%tuneResults = [1800 1160 1800 1800 1480 1160 1320 1800 1160 1480 1160 1480 1800 1640 1480 1480 1160 1640];
%CRFGM_test(leaveOut,iter,res);
end

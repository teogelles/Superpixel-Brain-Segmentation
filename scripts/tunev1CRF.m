function tune5CRF(leaveOut,res)
k = leaveOut; %silly way to avoid testing on training
while k == leaveOut
    k = randi(18,1); %hard coded for now
end
iter = CRFGM_v1Folds(k,1,1000,leaveOut,res);
fprintf('tuned to %d',iter);
%tuneResults = [1800 1160 1800 1800 1480 1160 1320 1800 1160 1480 1160 1480 1800 1640 1480 1480 1160 1640];
%CRFGM_test(leaveOut,iter,res);
end

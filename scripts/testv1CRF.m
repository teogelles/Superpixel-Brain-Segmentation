function testv1CRF(leaveOut,res)
k = leaveOut; %silly way to avoid testing on training
while k == leaveOut
    k = randi(18,1); %hard coded for now
end
starts = [600 600 600 600 800 800 800 800 600 600 600 600 900 900 900 500 500 500];
iter = CRFGM_v1test(k,1,starts(leaveOut),leaveOut,res);
fprintf('tuned to %d',iter);
%tuneResults = [1800 1160 1800 1800 1480 1160 1320 1800 1160 1480 1160 1480 1800 1640 1480 1480 1160 1640];
%CRFGM_test(leaveOut,iter,res);
end

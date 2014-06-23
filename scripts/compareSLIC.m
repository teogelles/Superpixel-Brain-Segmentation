function compareSLIC(res, numSuperVoxels, shapeParam, numIters)
    
    AD1Features = runSLIC(1, 500, 20, 15, 'AD', 1);
    AD2Features = runSLIC(1, 500, 20, 15, 'AD', 2);
    CN1Features = runSLIC(1, 500, 20, 15, 'CN', 1);
    CN2Features = runSLIC(1, 500, 20, 15, 'CN', 2);
    
    fprintf('AD1 Features:\n');
    disp(AD1Features);
    fprintf('\n');
    fprintf('AD2 Features:\n');
    disp(AD2Features);
    fprintf('\n');
    fprintf('CN1 Features:\n');
    disp(CN1Features);
    fprintf('\n');
    fprintf('CN2 Features:\n');
    disp(CN2Features);
    fprintf('\n');
    
end
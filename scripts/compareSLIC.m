function compareSLIC(res, numSuperVoxels, shapeParam, numIters)
    
    AD1Features = runSLIC(1, 'AD');
    AD2Features = runSLIC(2, 'AD');
    CN1Features = runSLIC(1, 'CN');
    CN2Features = runSLIC(2, 'CN');
    
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
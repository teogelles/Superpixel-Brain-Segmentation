function testEntropy(fileType, fileNum)
    
    numSVStart = 100;
    numSVStep = 20;
    numSVEnd = 500;

    SPStart = 0.05;
    SPStep = 0.05;
    SPEnd = 1;
    SPWeight = 20;
    
    entropyMatrix = zeros((numSVEnd-numSVStart)/numSVStep+1, ...
                           (SPEnd-SPStart)/SPStep+1);
    
    for i=numSVStart:numSVStep:numSVEnd
        for j=SPStart:SPStep:SPEnd
            
            avg = runSLIC(fileNum, fileType, 1, i, j);
            entropyMatrix(i, j*SPWeight) = avg;
        end
    end
    
    
    saveFileName = strcat('/scratch/tgelles1/summer2014/ADNI_Entropy/', ...
                          fileType, sprintf('%03d', fileNum))
    
    save(saveFileName, 'entropyMatrix')
end
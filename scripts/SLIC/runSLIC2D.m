function runSLIC2D(imageName, numSuperpixels, shapeParam, numIters) 
    
    imageAddr = strcat('/scratch/tgelles1/summer2014/crs/', ...
                       imageName, '.png');
    
    imageMatrix = imread(imageAddr);
    
    [labels border centerInfo] = SLIC_2D(imageMatrix,numSuperpixels, ...
                                         shapeParam, numIters);
    
    imwrite(uint8(labels), '/scratch/tgelles1/summer2014/temp/testLabels.png')
    imwrite(uint8(border), '/scratch/tgelles1/summer2014/temp/testBorders.png');
end
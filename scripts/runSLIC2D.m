function runSLIC2D(imageName, numSuperpixels, shapeParam, numIters) 
    
    imageAddr = strcat('/scratch/tgelles1/summer2014/crs/', ...
                       imageName, '.png');
    
    imageMatrix = imread(imageAddr);
    
    [labels border centerInfo] = SLIC_2D(imageMatrix,numSuperpixels, ...
                                         shapeParam, numIters);
    
    imwrite(labels, './temp/testLabels.png')
    imwrite(border, './temp/testBorders.png');
end
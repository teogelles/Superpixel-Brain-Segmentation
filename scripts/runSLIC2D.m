function runSLIC2D(imageAddr, numSuperpixels, shapeParam, numIters) 
    
    imageMatrix = imread(imageAddr);
    
    [labels border centerInfo] = SLIC_2D(imageMatrix,numSuperpixels, ...
                                         shapeParam, numIters);
    
    imwrite(labels, './temp/testLabels.png');
    imwrite(border, './temp/testBorder.png');
end
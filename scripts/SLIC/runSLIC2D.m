function runSLIC2D(imageName, numSuperpixels, shapeParam, numIters) 
    
    imageAddr = strcat('/scratch/tgelles1/summer2014/crs/', ...
                       imageName, '.png');

    fprintf('Loading image from: %s\n', imageAddr);
    imageMatrix = imread(imageAddr);
    
    [labels border centerInfo] = SLIC_2D(imageMatrix,numSuperpixels, ...
                                         shapeParam, numIters);
    
    labelAddr = '/scratch/tgelles1/summer2014/temp/testLabels.png';
    borderAddr = ['/scratch/tgelles1/summer2014/temp/' ...
                  'testBorders.png'];
    fprintf('Saving label to: %s\n', labelAddr);
    fprintf('Saving border to: %s\n', borderAddr);
    
    imwrite(uint8(labels), labelAddr);
    imwrite(uint8(border), borderAddr);
end
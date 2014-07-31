function testSCpipeSVM()
    
    accMatrix = zeros(6,6,10);
    
    fileName = ['/scratch/tgelles1/summer2014/slicExact120/' ...
                'features/CSV_NORM/organized_med.csv'];
    
    for i = 0:5
        for j = 0:5
            for k = 1:10
                clusters = 5+i*20;
                nebs = 5+j*20;
                sigma = k*.2;
                
                fprintf('%d clusters, %d neighbors, %2.1f sigma\n', ...
                        clusters,nebs,sigma);
                try
                    accMatrix(i,j,k) = SCpipeSVM(fileName,3,clusters,nebs,sigma);
                end
            end
        end
    end
    
    save('~/random/accMat.mat','accMatrix');
end
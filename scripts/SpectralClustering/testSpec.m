function testSpec()
    
    res = zeros(100,100,50);
    for i = 1:100
        for j = 1:100
            for k = 1:50
                try
                    SpectrallyCluster('organized_med.csv',i,j,k*.1)                
                    a = simpleMatch('med',1,2);
                    b = simpleMatch('med',2,2);
                    res(i,j,k) = (a+b)/2;
                end
            end
        end
    end
    save('~/random/res.mat',res);
    
    [maxVal, maxI] = max(res(:));
    
    [x y z] = ind2sub(size(res),maxI);
    
    fprintf('Best at %d clusters, %d neighbors, %2.1f sig\n',x,y,z* ...
            .1);
    fprintf('Best is %f\n',maxVal);
    
    figure
    surface(1:100,1:100,.1:5,res)
    
             
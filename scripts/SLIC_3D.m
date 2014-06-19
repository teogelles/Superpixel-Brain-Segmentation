%% SLIC, A Grayscale 3D MATLAB implementation
% The following code is a 3D MATLAB implementation of the SLIC
% algorithm, an elegantly simple algo used for image segmentation
% into a given number of superpixels following a given compatness
% parameter
%
% Currently it only works with black and white images (such as a
% matrix version of MR images)
%
% The original paper documenting this method can be found here
% http://infoscience.epfl.ch/record/177415/files/Superpixel_PAMI2011-2.pdf
%
% Paper written by Radhakrishna Achanta, Appu Shaji, Kevin Smith,
% Aurelien Lucchi, Pascal Fua, and Sabine Susstrunk
% Implemented by Andrew Gilchrist-Scott and Teo Gelles


function [labels,borders] = SLIC_3D(imageMat, numSuperVoxels, shapeParam)
    
    if ~(shapeParam) || (shapeParam < 0)
        shapeParam = 20;
        fprintf('Setting shapeParam to default of 20');
    end
    if ~(numSuperVoxels) || (numSuperVoxels < 0)
        numSuperVoxels = 200;
        fprintf('Setting numSuperVoxels to defualt of 200');
    end
    
    numVoxels = size(imageMat,1)*size(imageMat,2)*size(imageMat,3);
    % in the original cpp code, they add an odd .5 to S, but we're
    % not sure why, so we're going to leave that off for now
    step = (numVoxels/numSuperVoxels)^(1/3);
    
    % Initialize superpixel centers and adjust to neighbor point of
    % lowest gradient
    centers = getSeeds(imageMat,step);
    centers = adjustSeeds(imageMat, centers);
    
    % Initialize labels and distance matrices
    labels = -1*ones(size(imageMat));
    distances = Inf*ones(size(imageMat));
    
    for iterations = 1:10
        
        centerTracker = zeros(size(centers,1),5);
        
        for c = 1:size(centers,1)
            neb = getNeighborhood(imageMat,step,centers(c,1), ...
                                           centers(c,2),centers(c,3));
            for neb_i = 1:size(neb,1)
                D = calculateDistance(imageMat,centers(c,:), ...
                                      neb(neb_i,:),shapeParam, ...
                                      step);                
                if D < distances(neb(neb_i,1),neb(neb_i,2), ...
                                 neb(neb_i,3))
                    distances(neb(neb_i,1),neb(neb_i,2), ...
                              neb(neb_i,3)) = D;
                    labels(neb(neb_i,1),neb(neb_i,2), ...
                           neb(neb_i,3)) = c;
                    centerTracker(c,1) = centerTracker(c,1) + ...
                        neb(neb_i,1);
                    centerTracker(c,2) = centerTracker(c,2) + ...
                        neb(neb_i,2);
                    centerTracker(c,3) = centerTracker(c,3) + ...
                        neb(neb_i,3);
                    centerTracker(c,4) = centerTracker(c,4) + ...
                        imageMat(neb(neb_i,1),neb(neb_i,2),neb(neb_i,3));
                    centerTracker(c,5) = centerTracker(c,5) + 1;
                end
            end
        end
        
        newCenters = zeros(size(centers));
        
        parfor i = 1:size(centerTracker,1)
            if (centerTracker(i,5) == 0)
                newCenters(i,:) = centers(i,:);
                continue
            end
            
            newCenters(i,:) = centerTracker(i,1:4)./centerTracker(i, 5);
        end
        
        centers = newCenters;
        clear newCenters;
        
    end
    
    borders = getBorders(imageMat, labels, 1);

end

function seeds = getSeeds(imageMat, step)
% getSeeds takes the original image and the step size and returns a
% matriz of all the seed locations for starting the superpixel algo
    numSeeds = 0;
    n = 1;
    xstrips = int32(.5 + size(imageMat,1)/step);
    ystrips = int32(.5 + size(imageMat,2)/step);
    zstrips = int32(.5 + size(imageMat,3)/step);
    
    xerr = size(imageMat,1) - step*xstrips;
    if (xerr < 0)
        xstrips = xstrips - 1;
        xerr = size(imageMat,1) - step*xstrips;
    end
    
    yerr = size(imageMat,2) - step*ystrips;
    if (yerr < 0)
        ystrips = ystrips - 1;
        yerr = size(imageMat,2) - step*ystrips;
    end
    
    zerr = size(imageMat,3) - step*zstrips;
    if (zerr < 0)
        zstrips = zstrips - 1;
        zerr = size(imageMat,3) - step*zstrips;
    end
    
    xerrperstrip = xerr/xstrips;
    yerrperstrip = yerr/ystrips;
    zerrperstrip = zerr/zstrips;
    
    xoff = int32(step/2);
    yoff = int32(step/2);
    zoff = int32(step/2);
    
    numSeeds = xstrips*ystrips*zstrips;
    seeds = zeros(numSeeds,4);
    
    for z = 0:(zstrips-1)
        ze = z*zerrperstrip;
        d = z*step+zoff+ze;
        for y = 0:(ystrips-1)
            ye = y*yerrperstrip;
            for x = 0:(xstrips-1)
                xe = x*xerrperstrip;
                
                seeds(n,1:3) = [(x*step+xoff+xe) (y*step+yoff+ye) d];
                seeds(n,4) = imageMat(seeds(n,1),seeds(n,2),seeds(n,3));
                n = n+1;
            end
        end
    end
end

function seeds = adjustSeeds(im, seeds)
    grads = gradientApprox(im,seeds);
    for i = 1:size(seeds,1)
        ne = getNeighbors(im,seeds(i,1),seeds(i,2), seeds(i,3));
        for ne_i = 1:size(ne,1)
            if(abs(grads(ne(ne_i,1),ne(ne_i,2),ne(ne_i,3))) < ...
               abs(grads(seeds(i,1),seeds(i,2),seeds(i,3))))
                seeds(i,1:3) = ne(ne_i,:);
                seeds(i,4) = im(ne(ne_i,1),ne(ne_i,2),ne(ne_i,3));
            end
        end  
        clear ne
    end
    clear grads
end

function ne = getNeighbors(mat, i, j, k)
    num_ne = 1;
    % We claculate the number of neighbors so that we can
    % preallocate the space for the neighbors array
    if i == 1
        indi = [0 1];
        num_ne = num_ne * 2;
    elseif i == size(mat,1)
        indi = [-1 0];
        num_ne = num_ne * 2;
    else
        indi = [-1 0 1];
        num_ne = num_ne * 3;
    end
    if j == 1
        indj = [0 1];
        num_ne = num_ne * 2;
    elseif j == size(mat,2)
        indj = [-1 0];
        num_ne = num_ne * 2;
    else
        indj = [-1 0 1];
        num_ne = num_ne * 3;
    end
    if k == 1
        indk = [0 1];
        num_ne = num_ne * 2;
    elseif k == size(mat,3)
        indk = [-1 0];
        num_ne = num_ne * 2;
    else
        indk = [-1 0 1];
        num_ne = num_ne * 3;
    end
    
    ne = zeros(num_ne-1,3);
    ne_i = 1;
    
    for ii = indi
        for jj = indj
            for kk = indk
                if (ii == 0) && (jj == 0) && (kk == 0)
                    continue
                else
                    ne(ne_i,:) = [i+ii j+jj k+kk];
                    ne_i = ne_i + 1;
                end
            end
        end
    end              
end

function grads = gradientApprox(im,seeds)
%The following function approximates the gradient of the neighbors
%of the seeds values so that we can choose the lowest gradient
%value
    grads = zeros(size(im));
    
    for i = 1:size(seeds,1)
        
        if(grads(seeds(i,1),seeds(i,2),seeds(i,3)) ~= 0)
            continue
        end
        
        ne = getNeighbors(im,seeds(i,1),seeds(i,2),seeds(i,3));
        diffsum = 0;
        
        for ne_i = 1:size(ne,1)
            diffsum = diffsum + abs(im(seeds(i,1),seeds(i,2),seeds(i,3))...
                                    - im(ne(ne_i,1),ne(ne_i,2),ne(ne_i,3)));
            if(grads(ne(ne_i,1),ne(ne_i,2),ne(ne_i,3)) ~= 0)
                continue
            end
            
            nene = getNeighbors(im, ne(ne_i,1),ne(ne_i,2),ne(ne_i, ...
                                                             3));
            ne_diffsum = 0;
            
            for nene_i = 1:size(nene,1)
                ne_diffsum = ne_diffsum + ...
                    abs(im(ne(ne_i, 1),ne(ne_i,2),ne(ne_i,3)) - ...
                        im(nene(nene_i,1),nene(nene_i,2), ...
                           nene(nene_i,3)));
            end
            
            grads(ne(ne_i,1),ne(ne_i,2),ne(ne_i,3)) = ne_diffsum/ ...
                size(nene,1);
            clear nene
        end
        grads(seeds(i,1),seeds(i,2),seeds(i,3)) = diffsum/size(ne, 1);
        clear ne
    end
end

function ne = getNeighborhood(mat,sq_rad,i,j,k)
    num_ne = 1;
    % We claculate the number of neighbors so that we can
    % preallocate the space for the neighbors array
    if i-sq_rad <= 1
        indi = zeros(1,ceil(i+sq_rad));
        for ii = 1:(ceil(i+sq_rad))
            indi(ii) = ii;
        end
        num_ne = num_ne * ii;
    elseif i+sq_rad > size(mat,1)
        indi = zeros(1,size(mat,1) - floor(i-sq_rad));
        ind = 1;
        for ii = floor(i-sq_rad):size(mat,1)
            indi(ind) = ii;
            ind = ind+1;
        end
        num_ne = num_ne * (ind-1);
    else
        indi = floor(i - sq_rad):ceil(i+sq_rad);
        num_ne = num_ne * (2*ceil(sq_rad));
    end
    
    if j-sq_rad <= 1
        jndj = zeros(1,ceil(j+sq_rad));
        for jj = 1:(ceil(j+sq_rad))
            jndj(jj) = jj;
        end
        num_ne = num_ne * jj;
    elseif j+sq_rad > size(mat,2)
        jndj = zeros(1,size(mat,2) - floor(j-sq_rad));
        jnd = 1;
        for jj = floor(j-sq_rad):size(mat,2)
            jndj(jnd) = jj;
            jnd = jnd+1;
        end
        num_ne = num_ne * (jnd-1);
    else
        jndj = floor(j - sq_rad):ceil(j+sq_rad);
        num_ne = num_ne * (2*ceil(sq_rad));
    end
    
    if k-sq_rad <= 1
        kndk = zeros(1,ceil(k+sq_rad));
        for kk = 1:(ceil(k+sq_rad))
            kndk(kk) = kk;
        end
        num_ne = num_ne * kk;
    elseif k+sq_rad > size(mat,3)
        kndk = zeros(1,size(mat,3) - floor(k-sq_rad));
        knd = 1;
        for kk = floor(k-sq_rad):size(mat,3)
            kndk(knd) = kk;
            knd = knd+1;
        end
        num_ne = num_ne * (knd-1);
    else
        kndk = floor(k - sq_rad):ceil(k+sq_rad);
        num_ne = num_ne * (2*ceil(sq_rad));
    end
    
    ne = zeros(num_ne,3);
    ne_i = 1;
    
    for ii = indi
        for jj = jndj
            for kk = kndk
                % If the center is not wanted in the neighborhood,
                % uncomment this
                %if (ii == i) && (jj == j) && (kk == k)
                %    continue
                %end
                ne(ne_i,:) = [ii jj kk];
                ne_i = ne_i + 1;
            end
        end
    end
end

function dist = calculateDistance(mat,cent,neb,m,s)
    dsq = (cent(1)-neb(1))^2 + (cent(2)-neb(2))^2 + (cent(3)-neb(3))^2;
    dcq = (cent(4)-mat(neb(1),neb(2),neb(3)))^2;
    dist = sqrt(dcq + (dsq/(s^2))*(m^2));
end

function borders = getBorders(im,labels,fillSetter)
% If fill = 0, set the borders on the image to 0, else set them to
% inf
    if ~fillSetter
        fill = 0;
    else
        fill = inf;
    end
    
    borders = im;
    
    % We don't care about borders on the edge of the image, so we
    % start one voxel in
    for i = 2:(size(labels, 1)-1)
        for j = 2:(size(labels, 2)-1)
            for k = 2:(size(labels, 3)-1)
                breaker = 0;
                for sub_i = [i-1 i+1]
                    if labels(sub_i,j,k) ~= labels(i,j,k)
                        borders(i,j,k) = fill;
                        breaker = 1;
                        break
                    end
                end
                if breaker
                    continue
                end
                for sub_j = [j-1 j+1]
                    if labels(i,sub_j,k) ~= labels(i,j,k)
                        borders(i,j,k) = fill;
                        breaker = 1;
                        break
                    end
                end
                if breaker
                    continue
                end
                for sub_k = [k-1 k+1]
                    if labels(i,j,sub_k) ~= labels(i,j,k)
                        borders(i,j,k) = fill;
                        break
                    end
                end                                   
            end
        end
    end
end

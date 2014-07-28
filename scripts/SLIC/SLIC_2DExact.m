%% SLIC, A Grayscale 3D MATLAB implementation
% The following code is a 3D MATLAB implementation of the SLIC
% algorithm, an elegantly simple algo used for image segmentation
% into a given number of superpixels following a given compactness
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

function [labels, borders, centerInfo] = SLIC_2DExact(imageMat, ...
                                                 shapeParam, numIters)
% SLIC_2D - Get a superpixelated image
%
% @param imageMat - The image to be superpixelated as a matrix
% @param numSuperPixels - The number of superpixels to use (approximate)
% @param shapeParam - Weight used to change the importance of
% Euclidean distance in the calculateDistance function (which uses
% both Euclidean distance and a color metric)
% @param numIters - The number of iterations to run the SLIC loop
%
% @return lables - The labels for the superpixels of imageMat as a matrix
% @return borders - The boders for the superpixelated image overlayed on
% the original image
% @return centerInfo - Information on every superPixel including
% its average x,y, and z coordinates as well as the number of
% pixels that are members 
    
    if (shapeParam < 0)
        shapeParam = 20;
        fprintf('Setting shapeParam to default of 20');
    end

    imageMat = double(imageMat);
    
    numPixels = size(imageMat,1)*size(imageMat,2);

    % Initialize superpixel centers and adjust to neighbor point of
    % lowest gradient
    [centers steps] = getSeeds(imageMat);
    %centers = adjustSeeds(imageMat, centers);
    
    fprintf('Number of Centers Used: %d\n', size(centers, 1));
    
    % Initialize labels and distance matrices
    labels = -1*ones(size(imageMat));
    distances = Inf*ones(size(imageMat));
    
    
    imageMatSize = [size(imageMat,1),size(imageMat,2)];

    centerTracker = zeros(size(centers,1),4);

    normConst = max(imageMat(:));
    
    maxStep = max(steps);
    
    fprintf('Superpixelating Image');

    ijk = 1;
    
    totComps = numIters*size(centers, 1);
    tenPercent = totComps / 10;
    onePercent = totComps / 100;
    printCount = 1;
    for iterations = 1:numIters
        for c = 1:size(centers,1)
            
            if (printCount == tenPercent)
                fprintf('.');
                printCount = 1;
            else
                printCount = printCount + 1;
            end
            
            neb = getNeighborhoodEnds(imageMatSize,steps,centers(c,:));
            
            for i = neb(1):neb(2)
                for j = neb(3):neb(4)

                    curVox = [i j];
                    D = calculateDistance(imageMat,centers(c,:), ...
                                          curVox,shapeParam, ...
                                          maxStep, normConst);
                    

                    if ((D < distances(i,j)) || (labels(i,j) == c))
                        
                        distances(i, j) = D;
                        if  labels(i,j) ~= -1
                            centerTracker(labels(i,j),1) = ...
                                centerTracker(labels(i,j),1) ...
                                - i;
                            centerTracker(labels(i,j),2) = ...
                                centerTracker(labels(i,j),2) ...
                                - j;
                            centerTracker(labels(i,j),3) = ...
                                centerTracker(labels(i,j),3) ...
                                - imageMat(i,j);
                            centerTracker(labels(i,j),4) = ...
                                centerTracker(labels(i,j),4) ...
                                - 1;
                        end
                        
                        labels(i, j) = c;
                        centerTracker(c,1) = centerTracker(c,1) + i;
                        centerTracker(c,2) = centerTracker(c,2) + j;
                        centerTracker(c,3) = centerTracker(c,3) + ...
                            imageMat(i, j);
                        centerTracker(c,4) = centerTracker(c,4) + 1;
                    end
                end
            end
            ijk = 0;
        end
        
        % preallocation of the new centers
        newCenters = zeros(size(centers));
        for i = 1:size(centerTracker,1)
            if (centerTracker(i,4) == 0)
                newCenters(i,:) = centers(i,:);
                continue
            end
            
            newCenters(i,1) = centerTracker(i,1)/centerTracker(i, ...
                                                              4);
            newCenters(i,2) = centerTracker(i,2)/centerTracker(i, ...
                                                              4);
            newCenters(i,3) = centerTracker(i,3)/centerTracker(i, 4);
        end
        
        centers = newCenters;
    end

    fprintf('\n');
    
    borders = getBorders(imageMat, labels, 0);
    
    % build centerInfo
    centerInfo = zeros(size(centers,1),size(centers,2) + 1);
    centerInfo(:,1:size(centers,2)) = centers;
    centerInfo(:,size(centers,2) + 1) = centerTracker(:,end);
end

function [seeds steps] = getSeeds(imageMat)

    xDim = size(imageMat, 1);
    yDim = size(imageMat, 2);
        
    if (xDim < yDim)
        numSeeds = [21 24];
    else
        numSeeds = [24 21];
    end
    
    steps = [xDim/numSeeds(1) yDim/numSeeds(2)];
    
    seeds = zeros(numSeeds(1)*numSeeds(2), 3);
    
    index = 1;
    for i=1:numSeeds(1)
        for j=1:numSeeds(2)
            
            coords = round([(i-.5)*steps(1) (j-.5)*steps(2)]);
            
            neb = getNeighborhoodEnds([xDim yDim], [1 1], coords);
            
            neb = imageMat(neb(1):neb(2), neb(3):neb(4));
            
            avg = mean(neb(:));
            
            seeds(index, :) = [(i-.5)*steps(1) (j-.5)*steps(2) avg];
            
            index = index + 1;
        end
    end
end

function seeds = adjustSeeds(im, seeds)
% This function adjusts the seeds so that they are in the lowest
% gradient position of their immediate neighbors
% @param im - image matrix, shortened for ease of typing
% @param seeds - initial seeds
%
% @return seeds - seeds adjusted for gradient
    
    grads = gradientApprox(im,seeds);
    for i = 1:size(seeds,1)
        ne = getNeighbors(im,seeds(i,1),seeds(i,2));
        for ne_i = 1:size(ne,1)
            % if the gradient of the neighbor is less than that of
            % the current center, move the center there
            if(abs(grads(ne(ne_i,1),ne(ne_i,2))) < ...
               abs(grads(seeds(i,1),seeds(i,2))))
                seeds(i,1:2) = ne(ne_i,:);
                seeds(i,3) = im(ne(ne_i,1),ne(ne_i,2));
            end
        end  
        clear ne
    end
    clear grads
end

function ne = getNeighbors(mat, i, j)
% Gets the immediate neighbors of the i,j,k position in the matrix
% @param mat - image matrix
% @params i,j,k - x,y,z coordinates of the point of inters
%
% @return ne - number_of_neighbors x 3 matrix of neighbor
% coordinates
    
    num_ne = 1;
    % We calculate the number of neighbors so that we can
    % preallocate the space for the neighbors array

    % Edge case if we are on either edge, then else is anywhere else
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
    
    ne = zeros(num_ne-1,2);
    ne_i = 1;
    
    for ii = indi
        for jj = indj
            if (ii == 0) && (jj == 0)
                continue
            else
                ne(ne_i,:) = [i+ii j+jj];
                ne_i = ne_i + 1;
            end
        end
    end              
end

function grads = gradientApprox(im,seeds)
%The following function approximates the gradient of the neighbors
%of the seeds values so that we can choose the lowest gradient
%value
% @param im - image matrix
% @param seeds - initial seeded center coordinates
%
% @return grads = gradient approximation at each point that is a
% neighbor of a seed; the rest of the values we do not need to calculate
    
    grads = zeros(size(im));
    
    for i = 1:size(seeds,1)
        
        % if we already found the value, no need to do it again
        if(grads(seeds(i,1),seeds(i,2)) ~= 0)
            continue
        end
        
        % ne = neighbors
        ne = getNeighbors(im,seeds(i,1),seeds(i,2));
        diffsum = 0;
        
        for ne_i = 1:size(ne,1)
            % sum difference between neighboring value and center value
            diffsum = diffsum + abs(im(seeds(i,1),seeds(i,2))...
                                    - im(ne(ne_i,1),ne(ne_i,2)));

            if(grads(ne(ne_i,1),ne(ne_i,2)) ~= 0)
                continue
            end
            
            % nene = neighbors of neighbor
            nene = getNeighbors(im, ne(ne_i,1),ne(ne_i,2));
            
            ne_diffsum = 0;
            
            for nene_i = 1:size(nene,1)
                ne_diffsum = ne_diffsum + abs(im(ne(ne_i, 1),ne(ne_i,2)) ...
                                              - im(nene(nene_i,1), ...
                                                   nene(nene_i,2)));
            end
            
            % we average the difference by number of neighbors so
            % we don't get low approximations for edge cases
            grads(ne(ne_i,1),ne(ne_i,2)) = ne_diffsum/ ...
                size(nene,1);
            
            clear nene
        end
        grads(seeds(i,1),seeds(i,2)) = diffsum/size(ne, 1);
        clear ne
    end
end

function ne = getNeighborhood(mat,sq_rad,i,j)
    num_ne = 1;
    % We calculate the number of neighbors so that we can
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
    
    ne = zeros(num_ne,2);
    ne_i = 1;
    
    for ii = indi
        for jj = jndj
            ne(ne_i,:) = [ii jj];
            ne_i = ne_i + 1;
        end
    end
end

function dist = calculateDistance(mat,cent,neb,m,s, normConst)
% calculateDistance - Returns the distance between a pixel and its
% center with the special SLIC metric that incorporates both
% Euclidean distance and color
% @param mat - image matrix
% @param cent - center point, 1x4 matrix of x,y,z and intensity
% @param neb - neighbor to center, 1x3 matrix of x,y,z (intensity
% found in mat)
% @param m - the shape parameter we send in, generally in [1,40]
% @param s - the step size
    

    dsq = (cent(1)-neb(1))^2 + (cent(2)-neb(2))^2;
    % Square of Euclidean distance
        
    dcq = (cent(3)-mat(neb(1),neb(2)))^2;
    % Square of color distance
    
    distq = double(dcq/(normConst^2) + (dsq/(s^2))*(m^2));
    
    dist = sqrt(distq);
    % Overall distance
end



function neighborhoodEnds = getNeighborhoodEnds(imageMatSize, radius, coords)
% Function gets the neighborhood ends for a regions around the
% center of size radius
% @param imageMatSize - size of the image matrix
% @param radius - square radius of neighborhood around center
% @params i,j,k - x,y,z coordintate of center
%
% @return neighborhoodEnds - matrix of starts and ends of each
% direction of neighborhood
    
    i = coords(1);
    j = coords(2);
    
    neighborhoodEnds = [floor(i-radius(1)),ceil(i+radius(1)),floor(j-radius(2)), ...
                        ceil(j + radius(2))];
    
    % edge cases
    if neighborhoodEnds(1) < 1
        neighborhoodEnds(1) = 1;
    end
    if neighborhoodEnds(3) < 1
        neighborhoodEnds(3) = 1;
    end
    
    if neighborhoodEnds(2) > imageMatSize(1)
        neighborhoodEnds(2) = imageMatSize(1);
    end
    if neighborhoodEnds(4) > imageMatSize(2);
        neighborhoodEnds(4) = imageMatSize(2);
    end
end

function borders = getBorders(im,labels,fillSetter)
% Gets a matrix of the superpixel borders overlayed on the original
% image
% @param im - image matrix
% @param labels - final superpixel labels for each pixel
% @param fillsetter - If fill = 0, set the borders on the image to 0, else set them to
% inf
    if ~fillSetter
        fill = 0;
    else
        fill = inf;
    end
    
    borders = im;
    
    totIters = (size(labels, 1) - 2) * (size(labels, 2) - 2);
    
    iter10Percent = floor(totIters / 10);
    
    iterCounter = 0;
    % We don't care about borders on the edge of the image, so we
    % start one pixel in
    fprintf('Getting SLIC Border Overlay');
    for i = 2:(size(labels, 1)-1)
        for j = 2:(size(labels, 2)-1)
                
            iterCounter = iterCounter + 1;
            if (iterCounter == iter10Percent)
                fprintf('.')
                iterCounter = 0;
            end
            
            breaker = 0;
            for sub_i = [i-1 i+1]
                if labels(sub_i,j) ~= labels(i,j)
                    borders(i,j) = fill;
                    breaker = 1;
                    break
                end
            end
            if breaker
                continue
            end
            for sub_j = [j-1 j+1]
                if labels(i,sub_j) ~= labels(i,j)
                    borders(i,j) = fill;
                    breaker = 1;
                    break
                end
            end
            if breaker
                continue
            end
        end
    end
    
    fprintf('\n');
end
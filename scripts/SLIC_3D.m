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

function [labels, borders, centerInfo] = SLIC_3D(imageMat, numSuperVoxels, ...
                                                 shapeParam, numIters)
% SLIC_3D - Get a supervoxelated image
%
% @param imageMat - The image to be supervoxelated as a matrix
% @param numSuperVoxels - The number of supervoxels to use (approximate)
% @param shapeParam - Weight used to change the importance of
% Euclidean distance in the calculateDistance function (which uses
% both Euclidean distance and a color metric)
% @param numIters - The number of iterations to run the SLIC loop
%
% @return lables - The labels for the supervoxels of imageMat as a matrix
% @return borders - The borders for the supervoxelated image overlayed on
% the original image
% @return centerInfo - Information on every superVoxel including
% its average x,y, and z coordinates as well as the number of
% voxels that are members 
    
    if ~(shapeParam) || (shapeParam < 0)
        shapeParam = 20;
        fprintf('Setting shapeParam to default of 20');
    end
    if ~(numSuperVoxels) || (numSuperVoxels < 0)
        numSuperVoxels = 200;
        fprintf('Setting numSuperVoxels to defualt of 200');
    end
    
    numVoxels = size(imageMat,1)*size(imageMat,2)*size(imageMat,3);
    % in the original cpp code, they add an .5 to S. This is a
    % small trick that they do to get around rounding in int
    % conversion, however adding that to this code does slightly
    % change where the superVoxels are initialized and might be
    % worth applying. If you wish to do so, uncomment the line
    % below this and comment out the line below that
    %    step = .5 + (numVoxels/numSuperVoxels)^(1/3);
    step = round((numVoxels/numSuperVoxels)^(1/3));
    
    % Initialize superpixel centers and adjust to neighbor point of
    % lowest gradient
    centers = getSeeds(imageMat,step);
    centers = adjustSeeds(imageMat, centers);
    
    fprintf('Number of Centers Used: %d\n', size(centers, 1));
    
    % Initialize labels and distance matrices
    labels = -1*ones(size(imageMat));
    distances = Inf*ones(size(imageMat));
    
    
    imageMatSize = [size(imageMat,1),size(imageMat,2),size(imageMat,3)];
    
    fprintf('Supervoxelating Image');
    centerTracker = zeros(size(centers,1),5);
    % The algorithm technically calls for repeating this loop until
    % the change in placement of the centers is low, but as the
    % authors say 10 iterations generally suffices.  We use an
    % adjustable amount
    
    %disp(centerTracker);
    %disp(centers);
    for iterations = 1:2 %numIters
                         % disp(centers)
        
        fprintf('.');
        % centerTracker will keep track of the sum values in each
        % of the supervoxels so that at the end we can adjust the centers
        centerTracker(:,1:4) = 0;
        % centerTracker = zeros(size(centers,1),5);
        
        counter = 0;
        
        for c = 1:size(centers,1)
            
            neb = getNeighborhoodEnds(imageMatSize,step,centers(c,1), ...
                                                   centers(c, 2), ...
                                                   centers(c, 3));
            disp(neb)
            for i = neb(1):neb(2)
                for j = neb(3):neb(4)
                    for k = neb(5):neb(6)
                        counter = counter + 1;
                        curVox = [i j k];
                        D = calculateDistance(imageMat,centers(c,:), ...
                                              curVox,shapeParam, ...
                                              step);
                        if labels(i,j,k) == c
                            % disp('We are in!!!')
                        end
                        
                        if ((D < distances(i,j,k)) || ...
                                (labels(i,j,k) == c))
                            
                            distances(i, j, k) = D;
                            if  labels(i,j,k) ~= -1
                                centerTracker(labels(i,j,k),1) = ...
                                    centerTracker(labels(i,j,k),1) ...
                                    - i;
                                centerTracker(labels(i,j,k),2) = ...
                                    centerTracker(labels(i,j,k),2) ...
                                    - j;
                                centerTracker(labels(i,j,k),3) = ...
                                    centerTracker(labels(i,j,k),3) ...
                                    - k;
                                centerTracker(labels(i,j,k),4) = ...
                                    centerTracker(labels(i,j,k),4) ...
                                    - imageMat(i,j,k);
                                centerTracker(labels(i,j,k),5) = ...
                                    centerTracker(labels(i,j,k),5) ...
                                    - 1;
                            end
                            labels(i, j, k) = c;
                            %asdfsafsdfsa
                            % fprintf('in iter %d, setting center %d\n',iterations,c)
                            centerTracker(c,1) = centerTracker(c,1) + i;
                            centerTracker(c,2) = centerTracker(c,2) + j;
                            centerTracker(c,3) = centerTracker(c,3) + k;
                            centerTracker(c,4) = centerTracker(c,4) ...
                                + imageMat(i, j, k);
                            centerTracker(c,5) = centerTracker(c,5) + ...
                                1;
                        else
                            % disp(distances(i,j,k))
                        end
                    end
                end
            end
        end
        
        % disp(counter)
        
        % preallocation of the new centers
        newCenters = zeros(size(centers));
        
        % within this parfor loop (which can be changed back to
        % regular for if you don't have the appropriate MATLAB
        % toolbox) we assign the new centers to the indexes of the
        % old centers. If you wish to calculate the L2 norm of the
        % movement of the centers, E, you could do so here
        for i = 1:size(centerTracker,1)
            if (centerTracker(i,5) == 0)
                newCenters(i,:) = centers(i,:);
                continue
            end
            
            newCenters(i,1) = centerTracker(i,1)/centerTracker(i, ...
                                                              5);
            newCenters(i,2) = centerTracker(i,2)/centerTracker(i, ...
                                                              5);
            newCenters(i,3) = centerTracker(i,3)/centerTracker(i, ...
                                                              5);
            newCenters(i,4) = centerTracker(i,4)/centerTracker(i, 5);
        end
        
        centers = newCenters;
    end

    fprintf('\n');
    
    borders = getBorders(imageMat, labels, 0);
    
    % build centerInfo
    % disp(centerTracker);
    %disp(centers);
    centerInfo = zeros(size(centers,1),size(centers,2) + 1);
    centerInfo(:,1:size(centers,2)) = centers;
    centerInfo(:,size(centers,2) + 1) = centerTracker(:,end);
    centerInfo;
    %dsaknbad;kjbsalkhbvdkjsanbc;ljsaknfdlksajfdaslkfhsalf;khsafkhfsajfdhsadkfjsaf
    
end

function seeds = getSeeds(imageMat, step)
% getSeeds takes the original image and the step size and returns a
% matriz of all the seed locations for starting the superpixel algo
% @param imageMat = matrix of original 3D image
% @param step = step around each supervoxel that we want to compute
% our distances, this also encodes roughly our number of supervoxels
%
% @return seeds = the initial centers spread in a grid across the
% image matrix. Note that this will be approximately equal to the
% requested number of supervoxels, altered based upon the image dimensions
    
    numSeeds = 0;
    n = 1;
    % See remark in main code about the addition of .5
    %    xstrips = int32(.5 + size(imageMat,1)/step);
    %    ystrips = int32(.5 + size(imageMat,2)/step);
    %    zstrips = int32(.5 + size(imageMat,3)/step);    
    % number of superVoxels in each direction
    xstrips = round(size(imageMat,1)/step);
    ystrips = round(size(imageMat,2)/step);
    zstrips = round(size(imageMat,3)/step);
    
    % check that we don't have too many, if so adjust
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
    
    % find the number of voxels left off in each strip
    xerrperstrip = xerr/xstrips;
    yerrperstrip = yerr/ystrips;
    zerrperstrip = zerr/zstrips;
    
    % the initial offset in each direction
    xoff = int32(step/2);
    yoff = int32(step/2);
    zoff = int32(step/2);
    
    numSeeds = xstrips*ystrips*zstrips;
    seeds = zeros(numSeeds,4);
    
    % place the seeds
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
% This function adjusts the seeds so that they are in the lowest
% gradient position of their immediate neighbors
% @param im - image matrix, shortened for ease of typing
% @param seeds - initial seeds
%
% @return seeds - seeds adjusted for gradient
    
    grads = gradientApprox(im,seeds);
    for i = 1:size(seeds,1)
        ne = getNeighbors(im,seeds(i,1),seeds(i,2), seeds(i,3));
        for ne_i = 1:size(ne,1)
            % if the gradient of the neighbor is less than that of
            % the current center, move the center there
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
% @param im - image matrix
% @param seeds - initial seeded center coordinates
%
% @return grads = gradient approximation at each point that is a
% neighbor of a seed; the rest of the values we do not need to calculate
    
    grads = zeros(size(im));
    
    for i = 1:size(seeds,1)
        
        % if we already found the value, no need to do it again
        if(grads(seeds(i,1),seeds(i,2),seeds(i,3)) ~= 0)
            continue
        end
        
        % ne = neighbors
        ne = getNeighbors(im,seeds(i,1),seeds(i,2),seeds(i,3));
        diffsum = 0;
        
        for ne_i = 1:size(ne,1)
            % sum difference between neighboring value and center value
            diffsum = diffsum + abs(im(seeds(i,1),seeds(i,2),seeds(i,3))...
                                    - im(ne(ne_i,1),ne(ne_i,2),ne(ne_i,3)));

            if(grads(ne(ne_i,1),ne(ne_i,2),ne(ne_i,3)) ~= 0)
                continue
            end
            
            % nene = neighbors of neighbor
            nene = getNeighbors(im, ne(ne_i,1),ne(ne_i,2),ne(ne_i, ...
                                                             3));
            ne_diffsum = 0;
            
            for nene_i = 1:size(nene,1)
                ne_diffsum = ne_diffsum + ...
                    abs(im(ne(ne_i, 1),ne(ne_i,2),ne(ne_i,3)) - ...
                        im(nene(nene_i,1),nene(nene_i,2), ...
                           nene(nene_i,3)));
            end
            
            % we average the difference by number of neighbors so
            % we don't get low approximations for edge cases
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
% calculateDistance - Returns the distance between a pixel and its
% center with the special SLIC metric that incorporates both
% Euclidean distance and color
% @param mat - image matrix
% @param cent - center point, 1x4 matrix of x,y,z and intensity
% @param neb - neighbor to center, 1x3 matrix of x,y,z (intensity
% found in mat)
% @param m - the shape parameter we send in, generally in [1,40]
% @param s - the step size
    

    dsq = (cent(1)-neb(1))^2 + (cent(2)-neb(2))^2 + (cent(3)-neb(3))^2; ...
    % Square of Euclidean distance
    
    dcq = (cent(4)-mat(neb(1),neb(2),neb(3)))^2;
    % Square of color distance
    
    
    dist = sqrt(dcq + (dsq/(s^2))*(m^2));
    % Overall distance
end



function neighborhoodEnds = getNeighborhoodEnds(imageMatSize, radius, i, j, k)
    % Function gets the neighborhood ends for a regions around the
    % center of size radius
    % @param imageMatSize - size of the image matrix
    % @param radius - square radius of neighborhood around center
    % @params i,j,k - x,y,z coordintate of center
    %
    % @return neighborhoodEnds - matrix of starts and ends of each
    % direction of neighborhood
    
    neighborhoodEnds = [floor(i-radius),ceil(i+radius),floor(j-radius), ...
                        ceil(j + radius), floor(k - radius), ceil(k + radius)];
    
    % edge cases
    if neighborhoodEnds(1) < 1
        neighborhoodEnds(1) = 1;
    end
    if neighborhoodEnds(3) < 1
        neighborhoodEnds(3) = 1;
    end
    if neighborhoodEnds(5) < 1
        neighborhoodEnds(5) = 1;
    end
    
    if neighborhoodEnds(2) > imageMatSize(1)
        neighborhoodEnds(2) = imageMatSize(1);
    end
    if neighborhoodEnds(4) > imageMatSize(2);
        neighborhoodEnds(4) = imageMatSize(2);
    end
    if neighborhoodEnds(6) > imageMatSize(3);
        neighborhoodEnds(6) = imageMatSize(3);
    end
    
    if neighborhoodEnds(4) < neighborhoodEnds(3)
        disp('BadBad')
        disp(imageMatSize)
        disp(radius)
        disp(i)
        disp(j)
        disp(k)
    end
    
end

function borders = getBorders(im,labels,fillSetter)
% Gets a matrix of the supervoxel borders overlayed on the original
% image
% @param im - image matrix
% @param labels - final supervoxel labels for each voxel
% @param fillsetter - If fill = 0, set the borders on the image to 0, else set them to
% inf
    if ~fillSetter
        fill = 0;
    else
        fill = inf;
    end
    
    borders = im;
    
    totIters = (size(labels, 1) - 2) * (size(labels, 2) - 2) * ...
        (size(labels, 3) - 2);
    
    iter10Percent = floor(totIters / 10);
    
    iterCounter = 0;
    % We don't care about borders on the edge of the image, so we
    % start one voxel in
    fprintf('Getting SLIC Border Overlay');
    for i = 2:(size(labels, 1)-1)
        for j = 2:(size(labels, 2)-1)
            for k = 2:(size(labels, 3)-1)
                
                iterCounter = iterCounter + 1;
                if (iterCounter == iter10Percent)
                    fprintf('.')
                    iterCounter = 0;
                end
                
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
    
    fprintf('\n');
end
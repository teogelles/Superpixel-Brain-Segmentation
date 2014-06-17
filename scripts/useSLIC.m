% Code for using SLIC and outputting an image if called for

function useSLIC(im,sz,m_compact,dim)

    im = single(im);
    im2 = im;

    seg = vl_slic(im,sz,m_compact);
    disp('size of im')
    disp(size(im))
    disp('size of seg')
    disp(size(seg))
    % so the problem here is that the implementation of SLIC that
    % vl_feat has will not do 3d superpixels, so we're just going
    % to have to rewrite the code or make our own implementation
    % where the first ten is the superpixel size and the other is the
    % compactness parameter
    if (dim == 2)
        b = zeros(size(seg,1),size(seg,2));
        % b stands for boundaries
        for i = 1:size(seg,1)
            for j = 1:size(seg,2)
                ne = getNeighbors2(seg, i, j);
                for ne_i = 1:size(ne,1)
                    if seg(i,j) ~= seg(ne(ne_i,1),ne(ne_i,2))
                        b(i,j) = 255;
                        im2(i,j) = 0;
                    end
                end
            end
        end
    else
        b = zeros(size(seg,1),size(seg,2),size(seg,3));
        for i = 1:size(seg,1)
            for j = 1:size(seg,2)
                for k = 1:size(seg,3)
                    ne = getNeighbors3(seg, i, j, k);
                    for ne_i = 1:size(ne,1)
                        if seg(i,j,k) ~= seg(ne(ne_i,1),ne(ne_i,2),...
                                             ne(ne_i,3))
                            b(i,j,k) = 255;
                            im2(i,j,k) = 0;
                        end
                    end
                end
            end
        end
    end
    
    image(b)
    colormap gray
    figure
    image(im2)
    colormap gray
    
end

function ne = getNeighbors2(mat, i, j)
    if i == 1
        if j == 1
            ne = [i j+1; i+1 j];
        elseif j == size(mat,2)
            ne = [i j-1; i+1 j];
        else
            ne = [i j+1; i j-1; i+1 j];
        end
    elseif i == size(mat, 1)
        if j == 1
            ne = [i j+1; i-1 j];
        elseif j == size(mat,2)
            ne = [i j-1; i-1 j];
        else
            ne = [i j+1; i j-1; i-1 j];
        end
    else
        if j == 1
            ne = [i j+1; i-1 j; i+1 j];
        elseif j == size(mat,2)
            ne = [i j-1; i-1 j; i+1 j];
        else
            ne = [i-1 j; i j-1; i j+1; i+1 j];
        end
    end
end

function ne = getNeighbors3(mat, i, j, k)
    num_ne = 0;
    if i == 1
        indi = [0 1];
        num_ne = num_ne + 2;
    elseif i == size(mat,1)
        indi = [-1 0];
        num_ne = num_ne + 2;
    else
        indi = [-1 0 1];
        num_ne = num_ne + 3;
    end
    if j == 1
        indj = [0 1];
        num_ne = num_ne + 2;
    elseif j == size(mat,2)
        indj = [-1 0];
        num_ne = num_ne + 2;
    else
        indj = [-1 0 1];
        num_ne = num_ne + 3;
    end
    if k == 1
        indk = [0 1];
        num_ne = num_ne + 2;
    elseif k == size(mat,3)
        indk = [-1 0];
        num_ne = num_ne + 2;
    else
        indk = [-1 0 1];
        num_ne = num_ne + 2;
    end
    
    ne = zeros(num_ne,3);
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


% function ne = getNeighbors(mat, i, j)
%     if i == 1
%         if j == 1
%             ne = [i j+1; i+1 j; i+1 j+1];
%         elseif j == size(mat,2)
%             ne = [i j-1; i+1 j; i+1 j-1];
%         else
%             ne = [i j+1; i j-1; i+1 j; i+1 j-1; i+1 j+1];
%         end
%     elseif i == size(mat, 1)
%         if j == 1
%             ne = [i j+1; i-1 j; i-1 j+1];
%         elseif j == size(mat,2)
%             ne = [i j-1; i-1 j; i-1 j-1];
%         else
%             ne = [i j+1; i j-1; i-1 j; i-1 j-1; i-1 j+1];
%         end
%     else
%         if j == 1
%             ne = [i j+1; i-1 j; i-1 j+1; i+1 j; i+1 j+1];
%         elseif j == size(mat,2)
%             ne = [i j-1; i-1 j; i-1 j-1; i+1 j; i+1 j-1];
%         else
%             ne = [i-1 j-1; i-1 j; i-1 j+1; i j-1; ...
%                   i j+1; i+1 j-1; i+1 j; i+1 j+1];
%         end
%     end
% end
% Test code to see if the VLFeat libreary is worth anything

function testvl(sz,m_compact)

    close all;

    im = imread('IBSR_01_slice110.png');
    im = single(im);
    im2 = im;

    seg = vl_slic(im,sz,m_compact);
    % where the first ten is the superpixel size and the other is the
    % compactness parameter

    b = zeros(size(seg,1),size(seg,2));
    % b stands for boundaries
    for i = 1:size(seg,1)
        for j = 1:size(seg,2)
            ne = getNeighbors(seg, i, j);
            for ne_i = 1:size(ne,1)
                if seg(i,j) ~= seg(ne(ne_i,1),ne(ne_i,2))
                    b(i,j) = 255;
                    im2(i,j) = 0;
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

function ne = getNeighbors(mat, i, j)
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

function ne = getNeighbors(mat, i, j, k)
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
    
    ne = zeros(num_ne,3)
    ne_i = 0
        
    for ii = indi
        for jj = indj
            for kk = indk
                if (ii == i) && (jj == j) && (kk == k)
                    continue
                else
                    ne(ne_i,:) = [ii jj kk];
                    ne_i = ne_1 + 1;
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
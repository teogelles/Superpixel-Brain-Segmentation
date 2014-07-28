%% cropBlack.m
% Authors: Teo Gelles
%          Andrew Gilchrist-Scott
%
% This file contains the single function cropBlack, which removes
% from a 3D matrix all extraneous 0's, which in the case of the
% matrix being an image represents a black buffer around the real
% image

function [croppedImage indexList] = cropBlack(imageMatrix)
% cropBlack - crop the black buffer out of a matrix
%
% @param imageMatrix - The 3D matrix to crop
%
% @returns croppedImage - @imageMatrix with the black cropped off

    for i = 1:size(imageMatrix, 1)
        
        if any(any(imageMatrix(i, :, :)))
            xStartIndex = i;
            break;
        end
    end
    
    for j = 1:size(imageMatrix, 2)
        if any(any(imageMatrix(:, j, :)))
            yStartIndex = j;
            break;
        end
    end
    
    for k = 1:size(imageMatrix, 3)
        if any(any(imageMatrix(:, :, k)))
            zStartIndex = k;
            break;
        end
    end

    
    
    for i = size(imageMatrix, 1):-1:1
        if any(any(imageMatrix(i, :, :)))
            xEndIndex = i;
            break;
        end
    end
    
    for j = size(imageMatrix, 2):-1:1
        if any(any(imageMatrix(:, j, :)))
            yEndIndex = j;
            break;
        end
    end
    
    for k = size(imageMatrix, 3):-1:1
        if any(any(imageMatrix(:, :, k)))
            zEndIndex = k;
            break;
        end
    end
    
    indexList = zeros(3, 2);
    indexList(1, 1) = xStartIndex;
    indexList(1, 2) = xEndIndex;
    indexList(2, 1) = yStartIndex;
    indexList(2, 2) = yEndIndex;
    indexList(3, 1) = zStartIndex;
    indexList(3, 2) = zEndIndex;
    
    croppedImage = imageMatrix(xStartIndex:xEndIndex, yStartIndex: ...
                               yEndIndex, zStartIndex:zEndIndex);
end
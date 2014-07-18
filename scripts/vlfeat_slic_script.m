function vlfeat_slic_script(imageName, regionSize, regularizer)
   
    imageAddr = strcat('/scratch/tgelles1/summer2014/crs/', ...
                       imageName, '.png');
    
    im = imread(imageAddr);
    
    segments = vl_slic(single(im), regionSize, regularizer);
    
    perim = true(size(im,1), size(im,2));
    for k = 1 : max(segments(:))
        regionK = segments == k;
        perimK = bwperim(regionK, 8);
        perim(perimK) = false;
    end

    %perim = uint8(cat(3,perim,perim,perim));
    
    perim = uint8(perim);
    
    finalImage = im .* perim;
    imshow(finalImage);
end
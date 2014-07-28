% simple program that adjusts images for our presentation

function adjustImage()
    direc = '/scratch/agilchr1/pics/ToAdjust/';
    savedirec = '/scratch/agilchr1/pics/Adjusted/';
    listing = dir(direc);
    
    for i = 1:length(listing)
        if listing(i).name == '.'
            continue
        end
        
        im = imread([direc listing(i).name]);
        im = im*3;
        imwrite(im,[savedirec listing(i).name],'png');
    end
end

        
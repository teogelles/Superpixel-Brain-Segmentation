function findMissing()
    direc = '/scratch/tgelles1/summer2014/slicExact120/features';
    
    listing = dir(direc);
    AD = zeros(92,1);
    MCI = zeros(203,1);
    CN = zeros(102,1);
    
    for i = 1:length(listing)
        if listing(i).name(1) == '.' || length(listing(i).name) < 9
            continue
        end
        if listing(i).name(1) == 'M'
            type = listing(i).name(1:3);
            num = str2num(listing(i).name(4:6));
        else
            type = listing(i).name(1:2);
            num = str2num(listing(i).name(3:5));
        end
        if strcmp(type,'MCI')
            MCI(num) = 1;
        elseif strcmp(type,'AD')
            AD(num) = 1;
        elseif strcmp(type,'CN')
            CN(num) = 1;
        end
    end
    
    disp('AD')
    find(~AD)
    disp('CN')
    find(~CN)
    disp('MCI')
    find(~MCI)
end

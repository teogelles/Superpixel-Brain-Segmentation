% makes a simple surface out of our entropy values
% must be run in real gui matlab, not -nodisplay

filebase = '/scratch/tgelles1/summer2014/ADNI_Entropy/';
listing = dir(filebase)

x = .05:.05:1;
y = 100:20:500;

entropy = zeros(length(y),length(x));

for i = 1:length(listing)
    if strcmp(listing(i).name(1),'.');
        continue
    end
    
    filename = strcat(filebase,listing(i).name);
    rawEntropy = load(filename);
    rawEntropy = rawEntropy.entropyMatrix;
    entInd = find(rawEntropy);
    flatEnt = rawEntropy(entInd);
    ent = reshape(flatEnt,length(y),length(x));
    entropy = entropy + ent;
end

entropy = entropy/i;

figure
surface(x,y,entropy);

    
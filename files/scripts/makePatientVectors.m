function makePatientVectors(direc,numSV,alterCols,doMCI,equalAmounts)
% This is a simple matlab script that makes a matrix that has every
% patient as it's own vector with all the meta informaiton and the
% feature information from every SV
% @param - direc = directory where all the ANDI files that are to
% be concatenated exist
% @param - numSV = number of supervoxels in the images
% @param - alterCols = 
    
    if ~exist('numSV','var')
        numSV = 125;
    end
    
    if alterCols
        cols = [7:10 12 14];
        numFeat = length(cols);
    else
        numFeat = 18;
    end
    
    if doMCI
        types = {'AD','MCI','CN'};
    else
        types = {'AD','CN'};
    end
    
    % overallocate then we'll shrink
    total = zeros(92+203+102,15+numSV*numFeat);
    IDs = zeros(92+203+102,1);
    
    col_i = 0;
    for type_i = 1:length(types)
        if equalAmounts
            numEnd = 90;
        else
            numEnd = 205;
        end
        
        for num = 1:numEnd
            % briefly changing this to 90 so we get equal parts of
            % all data
            % We know that there are fewer patients than this, but
            % we don't really care since we check all files to see
            % if they exist
                        
            csvfile = strcat(direc,'CSV/',types{type_i}, ...
                             sprintf('%03d',num), '.csv');
            txtfile = strcat(direc,'TXT/',types{type_i}, sprintf('%03d', ...
                                                              num),'.txt');
            
            if ~exist(txtfile,'file')
                continue
            end
            
            row_i = 0;
            col_i = col_i + 1;            
            %input type
            IDs(col_i) = type_i;
            
            txtid = fopen(txtfile);
            txt = textscan(txtid,'%s','delimiter','\n');
            
            for i = 2:16
                
                lin = txt{1}{i,1};
                %This will capture all the metadata
                comInd = find(lin == ',');
                parenInd = find(lin == ')');
                
                numOfInterest = str2num(lin((comInd+1):(parenInd-1)));
                
                row_i = row_i + 1;
                
                total(col_i,row_i) = numOfInterest;
            end
            
            % Acts as basically a stus bar telling us how quickly
            % we're processing files
            disp(csvfile)
            csv = csvread(csvfile);
            
            %remove xyz
            if alterCols
                csv = csv(:,cols);
            end
            
            for i = 1:size(csv,1)
                total(col_i,(row_i+1):(row_i+numFeat)) = csv(i,:);
                row_i = row_i + numFeat;
            end
        end
    end

    total = total(1:col_i,:);
    IDs = IDs(1:col_i,:);
    
    total = normalizeData(total);
    
    savefilename = direc;
    savegroupname = direc;
    
    if equalAmounts
        savefilename = strcat(savefilename,'Eq');
        savegroupname = strcat(savegroupname,'Eq');
    end
    
    if alterCols
        savefilename = strcat(savefilename,'Intensity');
        savegroupname = strcat(savegroupname,'Intensity');
    end
    
    if doMCI
        savefilename = strcat(savefilename,'AllPat');
        savegroupname = strcat(savegroupname,'AllPat');
    else
        savefilename = strcat(savefilename,'ADCN');
        savegroupname = strcat(savegroupname,'ADCN');
    end
    
    fprintf('\nSaving files in:\n');
    savefilename = strcat(savefilename,'.csv')
    savegroupname = strcat(savegroupname,'_groups.csv')
    
    csvwrite(savefilename,total);
    csvwrite(savegroupname,IDs);
end

function total = normalizeData(total)
   
        total = (total - repmat(min(total,[],1), ...
                           size(total,1),1))*spdiags(1./ ...
                                                      (max(total, ...
                                                          [],1)-min(total,[],1))',0,size(total,2),size(total,2));
        goodCols = all(~isnan(total));
        total = total(:,goodCols);
end
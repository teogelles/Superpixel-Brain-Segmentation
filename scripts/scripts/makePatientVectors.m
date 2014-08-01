function makePatientVectors(direc,numSV)
% This is a simple matlab script that makes a matrix that has every
% patient as it's own vector with all the meta informaiton and the
% feature information from every SV
% @param - direc = directory where all the ANDI files that are to
% be concatenated exist
    
    if ~exist('numSV','var')
        numSV = 125;
    end
    
    alterCols = true;
    if alterCols
        cols = [7:10 12 14];
        numFeat = length(cols);
    else
        numFeat = 18;
    end
    
    types = {'AD','MCI','CN'};
    
    % overallocate then we'll shrink
    total = zeros(92+203+102,15+numSV*numFeat);
    IDs = zeros(92+203+102,1);
    
    col_i = 0;
    for type_i = 1:length(types)
        for num = 1:205
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
    
    savefilename = strcat(direc,'IntensityAllPat.csv');
    savegroupname = strcat(direc,'IntensityAllPat_groups.csv');
    csvwrite(savefilename,total);
    csvwrite(savegroupname,IDs);
end
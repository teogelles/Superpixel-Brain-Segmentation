% Testing and trying to develop an approximation algo for
% determining k factors of n which have the minimum spread between
% the min and the max factors

function testKFact(iStart,iEnd,kStart,kEnd)
    
    if ~exist('iEnd','var')
        iStart = 100;
        iEnd = 1000;
    end
    
    if ~exist('kEnd','var')
        kStart = 3;
        kEnd = 8;
    end
    
    types = {'minMin','minMed','minMax','buckets'};
    scores = [0, 0, 0, 0];
    times = [0, 0, 0, 0];

    for i = iStart:iEnd
        %fprintf('Num: %d\n',i);
        for j = 3:11
            tic;
            a = minMin(i,j);
            atime = toc;
            tic;
            b = minMed(i,j);
            btime = toc;
            tic;
            c = minMax(i,j);
            ctime = toc;
            tic;
            d = buckets(i,j);
            dtime = toc;
            
            times = times + [atime btime ctime btime];
            
            if (i == 504 || i == 336) && (j == 3);
                fprintf('For %d broken into %d:\n',i,j);
                disp(a);
                disp(b);
                disp(c);
                disp(d);
                fprintf('\n');
            end
            
            %check if algos worked at all
            if (prod(a) ~= i) || (length(a) ~= j)
                fprintf(['Problem with a; %d not %d prod, %d not %d ' ...
                         'factors\n'],prod(a),i,length(a),j);
            end
            if (prod(b) ~= i) || (length(b) ~= j)
                fprintf(['Problem with b; %d not %d prod, %d not %d ' ...
                         'factors\n'],prod(b),i,length(b),j);
            end
            if (prod(c) ~= i) || (length(c) ~= j)
                fprintf(['Problem with c; %d not %d prod, %d not %d ' ...
                         'factors\n'],prod(c),i,length(c),j);
            end 
            if (prod(d) ~= i) || (length(d) ~= j)
                fprintf(['Problem with d; %d not %d prod, %d not %d ' ...
                         'factors\n'],prod(d),i,length(d),j);
            end 
            
            spreadA = a(end) - a(1);
            spreadB = b(end) - b(1);
            spreadC = c(end) - c(1);
            spreadD = d(end) - d(1);
            spreads = [spreadA spreadB spreadC spreadD];
            
            if ~all(spreads == spreads(1))
                [~, ind] = min(spreads);
                type = types{ind};
                %fprintf('Winner: %s for %d broken in %d\n',type,i, j);
                scores(ind) = scores(ind) + 1;
            end
        end
    end
    
    times = times/(i*j);
    
    fprintf('Here are the scores:\n');
    for i = 1:length(scores)
        fprintf('%s : %d, avg time %.10f\n',types{i},scores(i),times(i));
    end
    
end

function k_facts = minMin(n,k)
    
    facts = factor(n);
    
    if length(facts) == k
        k_facts = facts;
        return
    elseif length(facts) < k
        k_facts = zeros(k,1);
        k_facts(1:(k-length(facts))) = 1;
        k_facts((k-length(facts)+1):end) = facts;
        return
    end
           
    while (length(facts) > k)
        newFact = facts(1)*facts(2);
        facts(2) = newFact;
        facts = sort(facts(2:end));
    end
    
    k_facts = facts;
end

function k_facts = minMax(n,k)
    
    facts = factor(n);
    
    if length(facts) == k
        k_facts = facts;
        return
    elseif length(facts) < k
        k_facts = zeros(k,1);
        k_facts(1:(k-length(facts))) = 1;
        k_facts((k-length(facts)+1):end) = facts;
        return
    end
    
    k_facts = [];
    while (length(k_facts) < k) && (length(facts) > 1) && ...
            (length([facts k_facts]) ~= k)
        newFact = facts(1)*facts(end);
        k_facts(end+1) = newFact;
        facts = facts(2:end-1);
    end
    
    if length([k_facts facts]) == k
        k_facts = [k_facts facts];
        k_facts = sort(k_facts);
    elseif length(k_facts) == k
        while(length(facts) ~= 0)
            k_facts = sort(k_facts);
            k_facts(1) = k_facts(1)*facts(1);
            facts = facts(2:end);
        end
        k_facts = sort(k_facts);
    else
        newFact = prod(facts);
        k_facts(end+1) = newFact;
        k_facts = sort(k_facts);
    end
end

function k_facts = minMed(n,k)
    
    facts = factor(n);
     
    if length(facts) == k
        k_facts = facts;
        return
    elseif length(facts) < k
        k_facts = zeros(k,1);
        k_facts(1:(k-length(facts))) = 1;
        k_facts((k-length(facts)+1):end) = facts;
        return
    end
        
    k_facts = [];
    while (length(k_facts) < k) && (length(facts) > 2) && ...
            (length([facts k_facts]) ~= k)
        newFact = facts(1)*facts(ceil(end/2));
        k_facts(end + 1) = newFact;
        facts = [facts(2:ceil(end/2)-1) facts(ceil(end/2)+1:end)];
    end
    
    if length([k_facts facts]) == k
        k_facts = [k_facts facts];
        k_facts = sort(k_facts);
    elseif length(k_facts) == k
        while(length(facts) ~= 0)
            k_facts = sort(k_facts);
            k_facts(1) = k_facts(1)*facts(1);
            facts = facts(2:end);
        end
        k_facts = sort(k_facts);
    else
        newFact = prod(facts);
        k_facts(end+1) = newFact;
        k_facts = sort(k_facts);
    end
    
end

function bucks = buckets(n,k)
    facts = factor(n);
    
    if length(facts) == k
        bucks = facts;
        return
    elseif length(facts) < k
        bucks = zeros(k,1);
        bucks(1:(k-length(facts))) = 1;
        bucks((k-length(facts)+1):end) = facts;
        return
    end
    
    bucks = ones(1,k);
    ideal = n^(1/k);
    
    randomize = false;
    if randomize
        randfactsind = randperm(length(facts));
        facts = facts(randfactsind);
    end
    
    reverse = false;
    if reverse
        facts = sort(facts,'descend');
    end
    
    for fact = facts
        err = zeros(size(bucks));
        for i = 1:length(bucks)
            err(i) = abs(ideal - bucks(i)*fact);
        end
        [~, minI] = min(err);
        bucks(minI) = bucks(minI)*fact;
    end
        
    % for fact = facts
    %     testVar = zeros(k,1);
    %     for i = 1:length(bucks)
    %         bucks(i) = bucks(i)*fact;
    %         testVar(i) = var(bucks);
    %         bucks(i) = bucks(i)/fact;
    %     end
    %     [~, minI] = min(testVar);
    %     bucks(minI) = bucks(minI)*fact;
    % end
        
    bucks = sort(bucks);
end
        
% Testing and trying to develop an approximation algo for
% determining k factors of n which have the minimum spread between
% the min and the max factors

function testKFact(iStart,iEnd)
    
    if ~exist('iEnd','var')
        iStart = 100;
        iEnd = 1000;
    end
    
    types = {'minMin','minMed','minMax'};
    scores = [0, 0, 0];

    for i = iStart:iEnd
        fprintf('Num: %d\n',i);
        for j = 3:8
            a = minMin(i,j);
            b = minMed(i,j);
            c = minMax(i,j);
            
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
            
            spreadA = max(a) - min(a);
            spreadB = max(b) - min(b);
            spreadC = max(a) - min(a);
            
            if (spreadA ~= spreadC) || (spreadB ~= spreadC) || ...
                    (spreadA ~= spreadB)
                [~, ind] = min([spreadA spreadB spreadC]);
                type = types{ind};
                fprintf('Winner: %s for %d broken in %d\n',type,i, j);
                scores(ind) = scores(ind) + 1;
            end
        end
    end
    
    scores
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
    while (length(k_facts) < k) && (length(facts) > 1)
        newFact = facts(1)*facts(end);
        k_facts(end+1) = newFact;
        fact = facts(2:end-1);
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
    
    k_facts = fact;
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
    while (length(k_facts) < k) && (length(facts) > 2)
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
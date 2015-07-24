%% Bootstrap
function aggs=bootstrap_revised(data,freqDist,Nsim)

freq = freqDist.random(Nsim, 1); 
m=memory;max_row=m.MaxPossibleArrayBytes/8.1;max_row=max_row/300;
% To be safe.  Dividing by a large nubmer will be safer.
freq_batch = cell(1); j = 0;  cf = uint64(cumsum(freq)); 

while  cf(end,1) > max_row;
    j = j + 1;    ending_index = find(cf < max_row,1,'last');
    freq_batch{j} = freq(1:ending_index);
    freq = freq((1+ending_index):end);   cf = uint64(cumsum(freq));
end
j = j + 1;    freq_batch{j} = freq;   aggs=[];


parfor k=1:length(freq_batch)  % Works even without Parallel Toolbox
 
     cf = cumsum(freq_batch{k});
     losses = randsample(data,cf(end),true);
     y = zeros(length(freq_batch{k}),1); 
     y(1,1) = sum(losses(1:(cf(1))),1);

    for j=2:length(cf)
    st=1+cf(j-1);en=cf(j);y(j,1)=sum(losses(st:en,1));
    end

   aggs = [aggs; y]; 
end    
 
end



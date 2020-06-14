%%%% Best Fit heuristic based approach of bin packing %%%%
function [z, sol] = BinpackFitness_BF(x, model)

n = model.n;
v = model.v;
Vmax = model.Vmax;
[~,idx] = sort(x);

%% Initialising empty bins and their max capacity
B = {};             
Bi = [];            
bin_cap = Vmax;     
rem_bin_cap = [];

%% Packing the items into the bins
Bi = [Bi,idx(1)];
B = [B;Bi];
bin_cap = bin_cap - v(idx(1));
rem_bin_cap = [rem_bin_cap;bin_cap];

for i=2:n
    % Identifying the bin with least remaining capacity that fits the
    % current object capacity
    [temp_val,temp_id] = sort(rem_bin_cap);
    Elg_bins = find(temp_val>v(idx(i)));
    if(numel(Elg_bins) ~= 0)
        %disp("Elgible bins: " + Elg_bins);
        j = temp_id(Elg_bins(1)); 
        Bi = B{j};
        Bi = [Bi,idx(i)];
        B{j} = Bi;
        rem_bin_cap(j) = rem_bin_cap(j) - v(idx(i));
    else  %% If none of the bins fits the current object a new bin is created and object is inserted into it
        Bi = [];
        Bi = [Bi,idx(i)];
        B = [B;Bi];
        bin_cap = Vmax - v(idx(i));
        rem_bin_cap = [rem_bin_cap;bin_cap];
    end
end

%% Computing the fitness of the solution obtained
nBin = numel(B);
Viol = zeros(nBin,1);
for i=1:nBin
    Vi = sum(v(B{i}));
    Viol(i) = Vi/Vmax;
end
fitness_cost = sum(Viol.^2) / nBin;

z = fitness_cost;
sol.nBin = nBin;
sol.B = B;
sol.percent_fill = Viol;
sol.fitnes = fitness_cost;
end
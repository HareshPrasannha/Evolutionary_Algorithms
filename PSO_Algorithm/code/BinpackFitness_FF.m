%%%% First Fit heuristic based approach of bin packing %%%%
function [z, sol] = BinpackFitness_FF(x, model)

n = model.n;
v = model.v;
Vmax = model.Vmax;

[~,idx] = sort(x,'descend');

%% Initialising empty bins and their max capacity
B = {};
Bi = [];
bin_cap = Vmax;

%% Packing items into the bin
for i=1:n
    % Identify if the current item fits the first available bin with free capacity
    if(v(idx(i)) < bin_cap)
        Bi = [Bi,idx(i)];
        bin_cap = bin_cap - v(idx(i));
	% If not create a new bin and pack the item
    else
        B = [B;Bi];
        %disp("Individual bins: " + Bi);
        Bi = [];
        bin_cap = Vmax;
        Bi = [Bi,idx(i)];
        bin_cap = bin_cap - v(idx(i));
    end
    if(i == n)
        B = [B;Bi];
    end
end

%% Computing the Fitness value of the solution obtained

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
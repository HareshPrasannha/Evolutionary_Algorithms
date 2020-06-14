%%%% Next Fit heuristic based approach of bin packing %%%%
function [z, sol] = BinpackFitness_NF(x, model)

n = model.n;
v = model.v;
Vmax = model.Vmax;

[~,idx] = sort(x);

%% Initialize empty bins and bin capacity
B = {};
Bi = [];
bin_cap = Vmax;
rem_bin_cap = [];

%% Packing items into the bin
Bi = [Bi,idx(1)];		%% Insert the first item into the bin
B = [B;Bi];
bin_cap = bin_cap - v(idx(1));	%% Update the remaining weight of the bin
rem_bin_cap = [rem_bin_cap;bin_cap];
flag = 1;

for i=2:n
	% Computing the next item that fits the available remaining bin capacity
    for j=1:numel(B)
        if(v(idx(i)) < rem_bin_cap(j))
            %disp("item number in inner for loop: "+idx(i));
            Bi = B{j};
            Bi = [Bi,idx(i)];
            B{j} = Bi;
            rem_bin_cap(j) = rem_bin_cap(j) - v(idx(i));
            flag = 0;
            break;
        end
    end
	% If none of the bins fit the current item create a new bin and insert the item
    if(flag == 1)
        %disp("Item number: " + idx(i));
        Bi = [];
        Bi = [Bi,idx(i)];
        B = [B;Bi];
        bin_cap = Vmax - v(idx(i));
        rem_bin_cap = [rem_bin_cap;bin_cap];
    end
    flag = 1;
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
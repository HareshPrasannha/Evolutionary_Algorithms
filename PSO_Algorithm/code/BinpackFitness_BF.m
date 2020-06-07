function [z, sol] = BinpackFitness_BF(x, model)
% clc;
% clearvars;
% close all;
% 
% model.v = [70,60,50,33,33,33,11,7,3];
% %[6 3 4 6 8 7 4 7 7 5 5 6 7 7 6 4 8 7 8 8 2 3 4 5 6 5 5 7 7 12];
% 
% model.n = numel(model.v);
% 
% model.Vmax = 100;
% 
% nVar = model.n;
% VarSize = [1 nVar];     % Decision Variables Matrix Size
% 
% VarMin = 0;     % Lower Bound of Decision Variables
% VarMax = 1;
% 
% x = unifrnd(VarMin,VarMax,VarSize);
%%
n = model.n;
v = model.v;
Vmax = model.Vmax;

[~,idx] = sort(x);

B = {};
Bi = [];
bin_cap = Vmax;
rem_bin_cap = [];
Bi = [Bi,idx(1)];
B = [B;Bi];
bin_cap = bin_cap - v(idx(1));
rem_bin_cap = [rem_bin_cap;bin_cap];
%%
for i=2:n
    %disp("item: " + v(idx(i)));
    %elg_bin_cap = 100;
    %elg_bin_id = 0;
    
    [temp_val,temp_id] = sort(rem_bin_cap);
    Elg_bins = find(temp_val>v(idx(i)));
    if(numel(Elg_bins) ~= 0)
        %disp("Elgible bins: " + Elg_bins);
        j = temp_id(Elg_bins(1));
        %j = elg_bin_id; 
        Bi = B{j};
        Bi = [Bi,idx(i)];
        B{j} = Bi;
        rem_bin_cap(j) = rem_bin_cap(j) - v(idx(i));
    else
        Bi = [];
        Bi = [Bi,idx(i)];
        B = [B;Bi];
        bin_cap = Vmax - v(idx(i));
        rem_bin_cap = [rem_bin_cap;bin_cap];
    end
    
%     for j=1:numel(B)
%         %disp("Bin"+j+": "+B{j});
%         if(v(idx(i)) < rem_bin_cap(j))
%             %disp("Eligible bin: " + j + "Rem_cap: " + rem_bin_cap(j));
%             if(elg_bin_cap > rem_bin_cap(j))
%                 elg_bin_cap  = rem_bin_cap(j);
%                 elg_bin_id = j;
%             end
%         end
%     end
%     if(elg_bin_id ~= 0)
%         %disp("Bin selected: " + elg_bin_id + "Selected bin cap: " + elg_bin_cap);
%         j = elg_bin_id; 
%         Bi = B{j};
%         Bi = [Bi,idx(i)];
%         B{j} = Bi;
%         rem_bin_cap(j) = rem_bin_cap(j) - v(idx(i));
%     else
%         %disp("Item number: " + idx(i));
%         Bi = [];
%         Bi = [Bi,idx(i)];
%         B = [B;Bi];
%         bin_cap = Vmax - v(idx(i));
%         rem_bin_cap = [rem_bin_cap;bin_cap];
%    end
end

nBin = numel(B);
Viol = zeros(nBin,1);

for i=1:nBin
    Vi = sum(v(B{i}));
    Viol(i) = Vi/Vmax;
end

fitness_cost = sum(Viol.^2) / nBin;
%cum_sum = 0;
% for i=1:numel(B)
%        cum_sum = cum_sum + (sum(v(B{i}))/ Vmax) ^ 2;
%        Fitness_regular = cum_sum/numel(B);
% end
z = fitness_cost;
sol.nBin = nBin;
sol.B = B;
sol.percent_fill = Viol;
sol.fitnes = fitness_cost;
%disp(sol.percent_fill);
end
function [z, sol] = BinpackFitness_NF(x, model)
% clc;
% clearvars;
% close all
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
flag = 1;
%%
for i=2:n
    for j=1:numel(B)
        %disp("Bin"+j+": "+B{j});
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
end
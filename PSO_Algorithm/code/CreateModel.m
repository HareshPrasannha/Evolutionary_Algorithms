function model = CreateModel()
    temp = csvread('u250_00.csv');     %% Set the path of the input dataset file appropriately
    model.v =  temp';                   %% Weights of different items to be packed in bins
    
    model.n = numel(model.v);           %% total no of items
    
    model.Vmax = 100;                   %% Max capacity of the bin

end
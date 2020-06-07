function model = CreateModel()
    temp = csvread('u250_00.csv');	%% Set the path of the input dataset file appropriately
    model.v =  temp';			%% Weights of different items to be packed in bins
    %[70,60,50,33,33,33,11,7,3];
    %[6 3 4 6 8 7 4 7 7 5 5 6 7 7 6 4 8 7 8 8 2 3 4 5 6 5 5 7 7 12];
    %[99,94,79,64,50,46,43,37,32,19,18,7,6,3];
    
    model.n = numel(model.v); 		%% total no of items
    
    model.Vmax = 150; 			%% Max capacity of the bin
    %100;
    %30;

end
%%%%% PSO FOR BIN PACKING ALGORITHM %%%%%
clc;
clear;  
close all;

%% Problem Definition

model = CreateModel();  % Bin Packing Problem initialisation (Set the path of the input dataset properly in the function "create model" appropiately

CostFunction = @(x) BinpackFitness_BF(x, model);  % Fitness Function

nVar = model.n;          % No of decision variables
VarSize = [1 nVar];      % Dimesion of Prticles

VarMin = 0;     % Lower Bound of Decision Variables
VarMax = 1;     % Upper Bound of Decision Variables

%% PSO Parameters

MaxIt=1000;      % No of Iterations

nPop=10;        % Swarm Size (No of particles in swarm)

% % PSO Parameters
% w=1;            % Inertia Weight
% wdamp=0.99;     % Inertia Weight Damping Ratio
% c1=1.5;         % Personal Learning Coefficient
% c2=2.0;         % Global Learning Coefficient

% Constriction Coefficients
phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt(phi^2-4*phi));
w=chi;          % Inertia Weight
wdamp=1;        % Inertia Weight Damping Ratio
c1=chi*phi1;    % Personal Learning Coefficient
c2=chi*phi2;    % Global Learning Coefficient

% Velocity Limits
VelMax=0.2*(VarMax-VarMin);
VelMin=-VelMax;

nParticleMutation = 8;       % Number of Mutations Performed on Each Particle
nGlobalBestMutation = 15;    % Number of Mutations Performed on Global Best

%% Initialization

empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Sol=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.Best.Sol=[];

particle=repmat(empty_particle,nPop,1);

GlobalBest.Cost=0;

for i=1:nPop
    
    % Initialize Position
    particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize);
    
    % Fitness Estimation
    [particle(i).Cost, particle(i).Sol]=CostFunction(particle(i).Position);
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    particle(i).Best.Sol=particle(i).Sol;
    
    % Update Global Best
    if particle(i).Best.Cost>GlobalBest.Cost
        GlobalBest=particle(i).Best;
    end
    
end

BestCost=zeros(MaxIt,1);
OptimalBins=zeros(MaxIt,1);

%% PSO Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        %IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        %particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Fitness Estimation
        [particle(i).Cost, particle(i).Sol] = CostFunction(particle(i).Position);
        
        % Perform Mutation
        for j=1:nParticleMutation
            NewParticle = particle(i);
            NewParticle.Position = Mutate(particle(i).Position);
            [NewParticle.Cost, NewParticle.Sol] = CostFunction(NewParticle.Position);
            if NewParticle.Cost >= particle(i).Cost
                particle(i) = NewParticle;
            end
        end
        
        % Update Personal Best
        if particle(i).Cost>particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            particle(i).Best.Sol=particle(i).Sol;
            
            % Update Global Best
            if particle(i).Best.Cost>GlobalBest.Cost
                GlobalBest=particle(i).Best;
            end
            
        end
        
    end
    
    % Perform Mutation on Global Best
    for i=1:nGlobalBestMutation
        NewParticle = GlobalBest;
        NewParticle.Position = Mutate(GlobalBest.Position);
        [NewParticle.Cost, NewParticle.Sol] = CostFunction(NewParticle.Position);
        if NewParticle.Cost >= GlobalBest.Cost
            GlobalBest = NewParticle;
        end
    end
    
    
    BestCost(it)=GlobalBest.Cost;       %% Storing global best fitness value
    OptimalBins(it)=GlobalBest.Sol.nBin;%% Storing global best bin pack solution
    
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it)) '  No of Bins = ' num2str(GlobalBest.Sol.nBin)]);
    
    w=w*wdamp;      %% Updating inertia weight factor
    
end

BestSol = GlobalBest;

writematrix(OptimalBins, "No.of_bins_Result.csv");
writematrix(BestCost, "Fitness_Result.csv");

%% Results

figure;
subplot(1,2,1),plot(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Fitness Value');
grid on;
subplot(1,2,2),plot(OptimalBins,'LineWidth',2);
xlabel('Iteration');
ylabel('Optimal No of Bins');
grid on;
suptitle("PSO for Bin packing problem Results");
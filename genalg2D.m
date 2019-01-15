function genalg2D

% 2-D Genetic Algorithm
%
% :: clear matlab's memory and close all open figures
clear, close all;

% :: some variables defining the population and its individuals
populationSize = 10; % the size of a population = number of individuals in a population
matrixSize = 10; % the size of the matrix we'll be using as a solution

% :: build a data structure to contain the initial population, which is a 3-dimensional array
population = zeros( matrixSize, matrixSize, populationSize );

% :: each individual in the population is initialized as a randomly filled matrixSize x matrixSize matrix
for( i = 1:populationSize )
    population( :, :, i ) = fix( 4 * rand( matrixSize, matrixSize ) );   
end

% :: a variable to keep track of the best fitness found so far
bestFitness = 0;

% :: calculate the optimal fitness
%     it is a node-edge problem so the maximum fitness can easily be computed
optimalFitness = 2 * ( ( matrixSize - 1 ) * matrixSize );

% :: let's keep track of the generation
generation = 0;

% :: do an infinite number of generations
while( 1 == 1 )
    
    % :: for each of the individuals in the population, find the fitness
    for( i = 1:populationSize )
        fitness( i ) = checkers( population( :, :, i ) );
    end
    
    % :: we're starting the new generation, keep track of this
    generation = generation + 1;
    
    % :: what's the best fitness in this population?    
    bestCurrentFitness = max( fitness );
    meanCurrentFitness = mean( fitness );
    
    % :: this will be used in generating a plot    
    plotFitnessMatrix( generation , 1 ) = generation;
    plotFitnessMatrix( generation , 2 ) = bestCurrentFitness;
    plotFitnessMatrix( generation , 3 ) = meanCurrentFitness;
    
    % :: if the optimal fitness is reached, break out    
    if( bestCurrentFitness == optimalFitness )
       disp(['Generation: ',num2str(generation)]);
       break
    end
    
    % :: displays every 25 generation output on the screen    
    if( mod(generation,25) == 0)
       disp(['Generation: ',num2str(generation),' ; Best Fitness:',num2str(bestCurrentFitness)]);
    end
        
    % :: following best fitness and plot the best solution so far on the screen
    if( bestCurrentFitness > bestFitness )
        % :: keep track of it        
        bestFitness = bestCurrentFitness;
        
        % :: show the best individual of the population in a figure
        %    this will only show the best individual if it's the best we've ever encountered
        showbestindividual( population, fitness );
    end
    
    % :: make some space to contain our newly build population    
    newPopulation = zeros( matrixSize, matrixSize, populationSize );
    
    % :: elitism: put the best inidvidual of the current population in the new population    
    newPopulation( :, :, 1 ) = bestindividual( population, fitness );
    
    % :: create children from the current population and put them in the new population    
    for( i = 2:populationSize )        
        % :: select an individual to act as parent        
        mom = select( population, fitness );
        
        % :: only use cross-over for 90% of the new children produced
        %    else copy the parent to the new population for 10% of the new children produced
        if( rand( 1 ) < 0.9 )         
            % :: select an extra individual to act as second parent            
            dad = select( population, fitness );
            
            % :: make a new child out of these parents            
            child = recombine( mutate( mom ), mutate( dad ) );
            
            % :: add the child to the new population            
            newPopulation( :, :, i ) = child;            
        else            
            newPopulation( :, :, i ) = mutate( mom );            
        end    
    end
    
    % :: overwrite the old population with the new one    
    population = newPopulation;
    
 end 

 % :: plot best individual to the screen 
 figure(1);
 showbestindividual( population, fitness );
 
 % :: plot the development of the fitness as a function of the time (generations)
 figure(2);
 plot(plotFitnessMatrix(:,1),[plotFitnessMatrix(:,2) plotFitnessMatrix(:,3)]);
 legend('Best Fitness','Mean Fitness',0);
 xlabel('Generations');
 ylabel('Fitness');
 
 
 
% --- THE SUBFUNCTION checkers - THIS IS AN INTERNAL FUNCTION ---
 
function f = checkers( individual )

% arguments: a matrix
% functionality: return the fitness of the coloring
[m, n] = size( individual );
f = 0;

for i = 1:m
    for j = 1:n
        k = individual( i, j );
        if( i < m )
            if( k ~= individual( i+1, j ) )
                f = f + 1;
            end
        end
        if( j < n )
           if( k ~= individual( i, j+1 ) )
               f = f + 1;
           end
        end
    end
end



% --- THE SUBFUNCTION showbestindividual - THIS IS AN INTERNAL FUNCTION ---

function showbestindividual( population, fitness )

% arguments: a population of individuals and a vector of fitnesses
% functionality: find the best individual in the population and show it in a figure
image( 25 * bestindividual( population, fitness ) ); % make a figure out of this best individual
colormap( jet ); % give it a color scheme
pause( .001 ); % give the computer some time to start drawing



% --- THE SUBFUNCTION bestindividual - - THIS IS AN INTERNAL FUNCTION ---

function result = bestindividual( population, fitness )

% arguments: a population of individuals and a vector of fitnesses
% functionality: find the best individual in the population
bestFitness = max( fitness ); % what is the best fitness?
bestIndividuals = find( fitness == bestFitness ); % which individuals have the best fitness
result = population( :, :, bestIndividuals(1) ); % we will only need one



% --- THE SUBFUNCTION select - - THIS IS AN INTERNAL FUNCTION ---

function selected = select( population, fitness )

% arguments: a population of individuals and a vector of fitnesses
% functionality: select an individual from the population, proportionate to the fitness
[m, n] = size( fitness );
s = sum( fitness );
roulette = s*rand( 1 );
currentSum = 0;
picked = n;

for i = 1:n
    if currentSum + fitness(i) < roulette
        picked = i;
    end
    currentSum = currentSum + fitness(i);
end

selected = population( :, :, picked );



% --- THE SUBFUNCTION recombine - - THIS IS AN INTERNAL FUNCTION ---

function child = recombine( mom, dad )

% arguments: two "parent" matrices
% functionality: return a recombination of these parents
sizeMom = size( mom );
sizeDad = size( dad );

if( sizeMom == sizeDad )
   
    child = zeros( sizeMom( 1 ) , sizeMom( 1 ) );
   
    % Selects a lower triangular matrix out of the Mom or Dad matrices and
    % combines them to a new child.
    
    diagParams = -sizeMom( 1 ) : 1 : sizeMom( 1 );
    sizeDiags = size( diagParams );
    indexDiagParam = fix( sizeDiags( 2 ) * ( rand( 1 ) ) ) + 1;
    cutDiagParam = diagParams( indexDiagParam );
   
    chooseMomOrDadAsUpper = fix( 2 * rand( 1 ) );
    
    if( chooseMomOrDadAsUpper == 0 ) % Mom as upper
       
       child = child + triu( mom , cutDiagParam );
       child = child + tril( dad , cutDiagParam - 1 );
              
    else % Dad as upper
       
       child = child + triu( dad , cutDiagParam );
       child = child + tril( mom , cutDiagParam - 1 );
      
    end
    
else % sizeMom != sizeDad
    
    disp( 'Error: when recombining, both parents should have the same size' );
    
end



% --- THE SUBFUNCTION mutate - - THIS IS AN INTERNAL FUNCTION ---

function mutated = mutate( original )

% arguments: an individual, a matrix
% functionality: return a mutated version of the matrix
[m n] = size( original );

for i = 1 : m
   for j = 1 : n
      if( rand( 1 ) > 0.99 )
         original( i , j ) = fix( 4 * rand( 1 ) );
      end
   end
end

mutated = original;   

% --- END OF CODE ---
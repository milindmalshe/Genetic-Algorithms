function GA1D(bitstringlength,population,mutation,cross_over,startleft);

% 1-D Genetic Algorithm
%
% Genetic algorithm to find the minimum of a function defined by the
% internal function FITFUNC.
%
% Input:
% - bitstringlength : the length of the string containing the bits
% - population      : the size of the population (needs to be divisible by two)
% - mutation        : probility that mutation happens (in percentage between 0 (0%) and 1 (100%))
% - cross_over      : probility that cross-over happens (in percentage between 0 (0%) and 1 (100%))
% - startleft       : Number of leading zeros in the bitstring. This overrules the random generation
%                      of bitstrings. This forces the population to be on the right side of the search space
%


% :: Checking input parameters
if ( nargin == 4 )
    startleft = 0;
elseif ( nargin < 4 )
    disp('Not enough input arguments');
else
    if ( startleft >= bitstringlength )
        disp('Too many leading zeros, reset to zero leading zeros');
        startleft = 0;
    else
        startleft = startleft;
    end
end
if ( ( mutation < 0 ) | ( mutation > 1 ) )
    disp('This mutationrate is not allowed. Reset the default value: 0.05 (5%)')
    mutation = 0.05;
elseif ( ( cross_over < 0 ) | cross_over > 1 )
    disp('This cross-overrate is not allowed. Reset the default value: 0.9 (90%)')
    cross_over = 0.9;
end
    
% :: Tracking numbers for the algorithm
k = 1;
bestfitness = max(fitfunc(0:1:(2^bitstringlength-1),bitstringlength));
meanfitness = [];
goalfitness = min(fitfunc(0:1:(2^bitstringlength-1),bitstringlength));

% :: Initialise the first population (random)
parents = (-1)*ones(bitstringlength+1,population);
parents(1:bitstringlength,1:population) = round(rand(bitstringlength,population));
parents(1:startleft,:) = 0;


% :: The MAIN loop -> The Genetic Algorithm
while ( bestfitness > goalfitness )

    % :: Determine fitness of the population and apend to parential matrix
    B = bin2dec(strcat(num2str(parents(1:bitstringlength,1:population))'));
    fitness = fitfunc(bin2dec(strcat(num2str(parents(1:bitstringlength,1:population))')),bitstringlength);
    parents(bitstringlength+1,:) = fitness';
    bestfitness = min(fitness);
    meanfitness = mean(fitness);

    % :: Sort the population according to its fitness
    [dump,index] = sort(parents(bitstringlength+1,:));
    sorted_parents = parents(1:bitstringlength,index);
    sorted_fitness = fitness(index);
    
    % :: Initialise new population
    new_population = (-1)*ones(bitstringlength,population);

    % :: Elitism (only and always only the best two)
    new_population(1:bitstringlength,1:2) = sorted_parents(:,1:2);
    
    % :: ROULETTE-WHEEL-SELECTION
    %     Select other members of new population and perform crossing-over and/or mutation
    relfit = ( 1 ./ sorted_fitness(1:1:population) );
    relfit = relfit / sum(relfit);
    
    % :: Going through reminder of the population to perform reproduction
    for i = 3:2:population

        % :: Determine which bitstring is the dad or the mom, based on the Roulette-Wheel
        dad_x = min(find(rand(1)<cumsum(relfit)));
        mom_x = min(find(rand(1)<cumsum(relfit)));
        
        % :: Select the two parents
        dad = sorted_parents(1:bitstringlength,dad_x);
        mom = sorted_parents(1:bitstringlength,mom_x);
    
        % :: Create childrens from parents
        child_1 = dad;
        child_2 = mom;
    
        % :: Perform crossing-over
        if ( rand(1) <= cross_over )
            cross_point = ceil(rand(1)*bitstringlength);
            child_1 = [dad(1:cross_point) ; mom(cross_point+1:bitstringlength) ];
            child_2 = [mom(1:cross_point) ; dad(cross_point+1:bitstringlength) ]; 
        end
    
        % :: Perform mutation
        children  = [ child_1 child_2 ];
        [row,col] = find(rand(bitstringlength,2)<=mutation);
        if ( isempty(row) == 0 )
            for j = 1:1:length(row)
                if ( children(row(j),col(j)) == 1 )
                    children(row(j),col(j)) = 0;
                else
                    children(row(j),col(j)) = 1;
                end
            end
        end
        child_1 = children(:,1);
        child_2 = children(:,2);
    
        % :: Put children in new population
        new_population(1:bitstringlength,i:i+1) = [ child_1 child_2 ];
        
    end

    % :: Make the new population the parents of the next generation
    parents = new_population;
    
    % :: Store fitness in an array and do the same with the generation
    b_f(k) = bestfitness;
    m_f(k) = meanfitness;
    gen(k) = k - 1;
    k = k + 1;
    
    % :: Generate some output to the screen
    figure(1);
    clf;
    plot(0:1:(2^bitstringlength-1),fitfunc(0:1:(2^bitstringlength-1),bitstringlength));
    v = axis;
    y_min = v(3);
    y_max = v(4);
    grid on;
    hold on;
    x = bin2dec(strcat(num2str(sorted_parents(1:bitstringlength,1:population))'));
    y = fitfunc(x,bitstringlength);
    for i = 2:1:population
        plot([x(i) x(i)],[y_min y_max],'g');
        plot(x(i),y(i),'xg');
    end
    plot([x(1) x(1)],[y_min y_max],'r');
    plot(x(1),y(1),'xr');
    text(20,0.3,sprintf('Generation %d',k-1));
    shg
    pause(.15);
    
end

% :: Plot the mean and the best fitness to the generations
figure(2);
plot(gen,[b_f ; m_f]);
xlabel('Generations');
ylabel('Fitness');
legend('Best Fitness','Mean Fitness');
grid on;
shg



% --- THE SUBFUNCTION FITFUNC - THIS IS AN INTERNAL FUNCTION ---

function y = fitfunc(x,bitstringlength)

factor = ( 2^bitstringlength - 1 );
y = cos(2*pi*x/factor).^2 + (4/5) * ( cos(pi*x/factor-pi/4) ) + (1/10) * ...
     ( sin(50*pi*x/factor - 25*pi) .* cos(10*pi*x/factor - 5*pi) ) + 1/2;
 
% --- END OF CODE ---
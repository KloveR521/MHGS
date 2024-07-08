
function [O1, O2]=UniformCrossover(i1,i2)

    alpha=randi([0 1],size(i1));
    
    O1=alpha.*i1+(1-alpha).*i2;
    O2=alpha.*i2+(1-alpha).*i1;
    
end
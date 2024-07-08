clear all %#ok<CLALL>
close all
clc

N=30; % Number of search agents   
T=1000; % Maximum number of iteration
D = 10;  
lb=-100;
ub=100;
MAIN_RUN = 30;
fhd=str2func('cec17_func');
for dd=1:30
    fprintf('\nCEC: %i\n', dd);
    F =strcat('CEC',num2str(dd));
    Function_name=F; % Name of the test function, range from F1-F13
    allrun_fitness=zeros(1,MAIN_RUN);
    func_num=dd;
    for mainloop=1:MAIN_RUN
        fprintf('\nCEC: %d Run: %i\n',dd,mainloop);
        
        [Destination_fitness,bestPositions,Convergence_curve]=RMHGS(fhd,D,N,T,lb,ub,func_num);

        allrun_fitness(mainloop)=Destination_fitness;
    end
    mean_fitness=mean(allrun_fitness);
    std_fitness=std(allrun_fitness);
    save(strcat('RMHGS_30_1000_F',num2str(dd)),'allrun_fitness','mean_fitness','std_fitness');
    
    result(dd).name=F;
    result(dd).mean_fitness=mean_fitness;
    result(dd).std_fitness= std_fitness;
    save('result_RMHGS_30_1000.mat','result');

end

%% C_ILD_standalone.m
% This models the regulatory capital using the selected ILD severity
% distributions.
%
clear all;
Internal_Data_A_load_data;

nSims = 1e6;   
Param = [];



for l = 1:length(UoMs)
    %% Load UoM specific parameters
    % Select the appropriate input parameters for this UoM
   
    um = l; 
    disp(['Testing UoM ' num2str(um)]);
    ILDThres = ILD_thresholds(l);
    if isnan(ILDThres) % use modeling threshold if no other threshold is specified
        ILD_thresholds(l) = modeling_threshold;
        ILDThres = ILD_thresholds(l);
    end
    ILD_dist = ILD_dists(l);
    
    % Load the ILD and ELD columns
    Port = D(:,1); Year = D(:,3); Qrt=D(:,4);   data_y = D(:,2); data_yb =D(:,2);
    
    % Set up function parameters
    ControlVars.GofTest.Sims = 500; %number of simulations for computing gof tests pvalue
    ControlVars.GofTest.Alpha = 0.05; %significance level for gof tests
    ControlVars.FitTruncatedDistributions_m.InitialFit.MaxIter = 5e3; %max iterations for optimization algorithm
    ControlVars.FitTruncatedDistributions_m.InitialFit.MaxFunEvals = 10e3;%max function evaluation for optimization algorithm
    ControlVars.FitTruncatedDistributions_m.GofFit.MaxIter = []; %max iterations for optimization algorithm used during gof tests (empty value means matlab default)
    ControlVars.FitTruncatedDistributions_m.GofFit.MaxFunEvals = [];%max function evaluation for optimization algorithm used during gof tests (empty value means matlab default)
    ControlVars.lThreshold.Internal_l = modeling_threshold;  % 2e4 tail threshold
    lThres_l = ControlVars.lThreshold.Internal_l;
    ControlVars.lThreshold.Internal = ILDThres;
    lThres = ControlVars.lThreshold.Internal;
    
    % Select the UoM data above threshold
    data = data_y(Port==um & data_y>=lThres);
    data_b = data_y(Port==um & data_y>=lThres_l & data_y<lThres);
    max_e=max(data);
    
 
    
    % Compute annual frequencies
    if um~=3
        
        event_y  = length(data_y(Port==um & data_y>=lThres & ((Year+Qrt/4)>=(year_max-5)) ))/5;
        event_yb = length(data_y(Port==um & data_y>=lThres_l & data_y<lThres & ((Year+Qrt/4)>=(year_max-5)) ))/5;
      
        
    else % CPBP frequencies are computed at 13 years
        event_y = length(data_y(Port==um & data_y>=lThres))/(year_max-year_start);
        event_yb = length(data_y(Port==um & data_y>=lThres_l & data_y<lThres ))/(year_max-year_start);
    end
    
    % Updated Lambda
%     if ~isnan(updated_lambda(l));
%         event_y_original  = event_y;
%         event_yb_original =  event_yb;
% 
%         event_y   = updated_lambda(l) * event_y_original  / (event_yb_original + event_y_original) ; 
%         event_yb  = updated_lambda(l) * event_yb_original / (event_yb_original + event_y_original) ; 
%     end
    
    
    if ~isnan(updated_lambda(l));
        if lThres == 2e4
            event_y = updated_lambda(l);
        else
           event_yb = updated_lambda(l) - event_y;
        end
    
    end
    
    
    
    freqDist_b = makedist('Poisson', 'lambda', event_yb);
    freqDist = makedist('Poisson', 'lambda', event_y);

    %% Bootstrap the body of the distribution 
    if lThres_l<lThres
        aggs_b=bootstrap_revised(data_b,freqDist_b,nSims);
    else
        aggs_b=0;
    end;
    
    boot_cap = quantile( aggs_b+ bootstrap_revised(data,freqDist,nSims),.999)/1e6;
    
    
    %% Produce the standalone ILD distributions
    distr=['Lognormal         ';'Weibull           ';'Loggamma          ';'Burr              ';'Loglogistic       ';'Generalized Pareto';'LgnormmixPSA      ';'LgnormmixRSA      '];
    cell_dis=cellstr(distr);
    disp('Fitting tail distribution...' )
    
    % Fit ILD distribution
   
    
    
    for i=1:8
        con_f=[];     disp(['Severity Distribution = ',num2str(i),' for UoM = ' num2str(um)]);
        fittedDist = cell(1,2);      distri=char(cell_dis(i));
        [ fittedDist{1}, fittedDist{2} ] = fitTruncatedDistributions_modified(data,distri,lThres,ControlVars);
      
        b_t=fittedDist{1,1}.cdf(lThres);   m=fittedDist{1,1}.mean;       v=fittedDist{1,1}.var;
               
        con_e= struct2cell(fittedDist{2});
        try;    con_f=con_e{1}.message;   catch;  con_f='0';     end;
                
        if  b_t < 1
           figure('visible','off')
           plot3_99_9percentiles(data, lThres, fittedDist{1});
           print('-djpeg','-r200',[pwd '\QQplots\QQ_Plot_UoM' num2str(um) '_Dist' num2str(i) '_thr' num2str(lThres) '.jpg']) 
        end
    
        %% Perform gof tests on the distribution
        disp('Performing GoF tests...')
        TestResult = gofTest_parallel_modified(fittedDist,'truncated',ControlVars);
        Test_sc=GoFscore(TestResult);
         
      disp('Performing Monte Carlo...')
    
    
    % Simulate losses above ILD threshold
   aggs = nan;  attloss = 0; %attritional loss (mean aggregate loss below the threshold)
    if  b_t < 1
      aggs = runMonteCarlo_revised(fittedDist{1}, lThres, freqDist,attloss, nSims);
    end
    
    disp(['1 in 1000 year loss is: ' num2str(round((quantile(aggs, 0.999)+quantile(aggs_b, 0.999))/1e6)) 'MM for UoM ' num2str(um)]);
    
    
        %% Output to excel
        Tmp=cellstr(char({ num2str(um), num2str(lThres),distri,num2str(TestResult.KS.pvalue), num2str(TestResult.AD.pvalue),...
            num2str(TestResult.AD2.pvalue), num2str(TestResult.CVM.pvalue), ...
            num2str(TestResult.ADup.pvalue), num2str(TestResult.AD2up.pvalue),num2str(Test_sc),...
            num2str(m),num2str(v),num2str(b_t),num2str(con_f),...
             num2str(max(data)),num2str(event_yb),num2str(event_y),...
             num2str((quantile(aggs, 0.999)+quantile(aggs_b, 0.999))/1e6), num2str(boot_cap),num2str(length(data))}));
        Param = [Param;Tmp'];
    end



end

%% Output to excel file
output_file = [output_folder  '\' datestr(now,'yyyymmdd_HHMM_') 'ILD_Dist_Fit_stats.xlsx'];
colnames = {'UOM','Threshold','Distribution', 'KS', 'AD', 'AD2', 'CVM', 'ADup', 'AD2up', 'Score',...
            'Mean', 'Variance', 'Below threshold','Convergence Problems', 'Max Loss',...
            'Lambda body','Lambda tail','Parametric Capital MM','Bootstrap Capital MM','Severity Tail Obs N'};

m = num2cell(Param);
outmatrix =cell2table([colnames; m]);
writetable(outmatrix,output_file,'WriteVariableNames',false);




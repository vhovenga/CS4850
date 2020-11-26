
function CBCR(input, curricula, alpha, gamma_1, gamma_2, learning_rate, max_iter, final_iter)
disp ('============================================');
disp ('Data = Hi-C data  ');
disp ('===========================================');
path = 'Scores/';
%==========================================================================
% Make directory if it doesn't exist
%--------------------------------------------------------------------------
if ~exist('Scores', 'dir')
    % Folder does not exist so create it.
    mkdir('Scores');
end
%--------------------------------------------------------------------------
% 3DMax Variables
%--------------------------------------------------------------------------

smooth_factor = 1e-6; % for numerical stability

alpha = alpha; % Scaling factor for trained data.
beta = 1-alpha; % Scaling factor for untrained data. 
gamma_1 = gamma_1; % Decay rate for first moment estimate in adam. 
gamma_2 = gamma_2; % Decay rate for second moment estimate in adam.
NUMBER_OF_CURRICULA = curricula;

NEAR_ZERO =0.00001; % used to signify a boundary of convergence

NUM = 1; 
  
LEARNING_RATE = learning_rate; % Specify the learning rate.

INPUT_FILE = input;
    

MAX_ITERATION = max_iter; % maximum number of iterations
FINAL_ITERATION = final_iter; % Maximum number of iteration on final training
%==========================================================================

[filepath, name, ext] = fileparts(INPUT_FILE);

S  = [];
TS = [];
Corr = [];
P_CORR = [];
RMSD = [];
new_name = name;


%=========================================================================
ReadInput; % Load the contact file
%=========================================================================



% specify the alpha


for CONVERT_FACTOR = .1:0.1:2


        %-------------------------
        Convert2Distance;
        %-------------------------
        
        %-------------------------
        Divide_Curricula; % Divide into curricula
        %-------------------------

        % For each CONVERT_FACTOR generate N structures
        s= [];
        ts=[];
        cor=[]; 
        P_corr = [];
        rmsd = [];
        totalIF = 0; 
        prev_str = [];
        l = 0;
        
        
        prev_trained_data = []; % Matrix all trained curricula 
        prev_str = []; % Vector of xyz data for all trained structures from each previous curriculum. 
        
        for l = 1:(1+NUMBER_OF_CURRICULA)
            %disp ('=============================================================================');
            %fprintf('Creating structure at CONVERT_FACTOR = %f and structure = %d\n', CONVERT_FACTOR,l);
            %disp ('=============================================================================');
            data = [C{l}; prev_trained_data];
            %-------------------------
            Optimization;
            %-------------------------
  
            prev_trained_data = data;
            
        end
            %========================================================================
            % scoring using spearman correlation, pearson correlation and  RMSD     
            %------------------------------------------------------------------------
             str_name =[ path,num2str(name),'_CONVERT_FACTOR=',num2str(CONVERT_FACTOR),'N=',num2str(l)];
             %----------------------
             Evaluation;
             %-----------------------             
             %output scores ::: S shows the scores for different CONVERT_FACTOR            
             ts = [ts;CONVERT_FACTOR,l,SpearmanRHO];                        
             cor = [cor;SpearmanRHO];  %output the correlation         
            %----------------------
            %pearson correlation           
            P_corr = [P_corr;PearsonRHO];
            %----------------------
            % perform root mean square error
            rmsd = [rmsd;rmse];
            %------------------------------------------------------------------------     
            %output pdb and image
            %------------------------------------------------------------------------
            output; % default run
  
         %------------------------------------------------------------------------     
         %Select the representative model for the CONVERT_FACTOR value
         %------------------------------------------------------------------------
           [v,index]= max(ts(:,3));
           TS = [TS;ts(index,:)];      
           Corr = [Corr;max(cor)];
           P_CORR = [P_CORR ;max(P_corr)];
           RMSD = [RMSD;min(rmsd)];    
end

% ========================================================================
% Select the representative model for the chromosome
%---=--========------------------------------------------------------------
[first, name, ext] = fileparts(INPUT_FILE);
input = strcat(name,ext);
[val,index]= max(TS(:,3));
CONVERT_FACTOR = TS(index,1);
l = TS(index,2);
pl =  P_CORR(index);
rep_rms = RMSD(index);
% save scores to file
f_scores =strcat(path,num2str(new_name),'_Finalscores.txt');
%dlmwrite(f_scores,sprintf('CONVERT_FACTOR\tStructure Number\tSpearman Correlation'));
%dlmwrite(f_scores,TS);
%save pearson correlation
scores =strcat(path,num2str(new_name),'_pearsoncorr.txt');
%dlmwrite(scores, P_CORR );
%save spearman correlation
scores =strcat(path,num2str(new_name),'_spearmancorr.txt');
%dlmwrite(scores,Corr);
%save RMSD
scores =strcat(path,num2str(new_name),'_rmsd.txt');
%dlmwrite(scores,RMSD);

final_name =[ num2str(new_name),'_CONVERT_FACTOR=',num2str(CONVERT_FACTOR),'N=','.pdb']; 
readme =strcat(path,num2str(new_name),'_readme.txt');
fprintf('\n===============================================\n');
fprintf('Input file: %s\n', input)
fprintf('Convert factor: %f\n', CONVERT_FACTOR);  
fprintf('AVG RMSE: %f\n', rep_rms)
fprintf('AVG Spearman correlation Dist vs. Reconstructed Dist: %f\n', val);  
fprintf('AVG Pearson correlation Dist vs. Reconstructed Dist: %f\n', pl);  
%fid = fopen(readme,'wt');
%msg = ['The asdfadfi structure is ',final_name]; 
%fprintf(fid,msg);
%fclose(fid);

out_log = strcat(path, num2str(new_name),'_output.log');
log = fopen(out_log,'wt');
fprintf(log,'Input file: %s\n', input);
fprintf(log, 'The optimal structure is: %s\n', final_name);
fprintf(log,'Convert factor: %f\n', CONVERT_FACTOR);  
fprintf(log,'AVG RMSE: %f\n', rep_rms);
fprintf(log,'AVG Spearman correlation Dist vs. Reconstructed Dist: %f\n', val);  
fprintf(log,'AVG Pearson correlation Dist vs. Reconstructed Dist: %f\n', pl); 
fclose(log);
end 


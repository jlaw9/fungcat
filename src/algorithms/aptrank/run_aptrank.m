function runAptRank
% Set up cvx
% (1) Download cvx from http://cvxr.com/cvx/download/ and unzip the package.
% (2) Add the cvx path
addpath('/home/jeffl/git-workspace/aptrank/cvx');
% (3) Set it up:
cvx_setup

% Load the datasets:
%network_file = '/data/jeff-law/projects/fungcat-function-prediction/test/2017_10-seq-sim-x5-string-net.mat';
%annotations_and_go_dag_file = '/data/jeff-law/projects/fungcat-function-prediction/test/P-annotations-and-go-dag.mat';
network_file = '/data/jeff-law/projects/fungcat-function-prediction/test/aptrank2/2017_10-string-208964-net.mat';
annotations_and_go_dag_file = '/data/jeff-law/projects/fungcat-function-prediction/test/aptrank2/P-annotations-and-go-dag.mat';
% contains G
load(network_file);
% contains R and H for a given GO category (i.e., biological process)
load(annotations_and_go_dag_file);
% for some reason all of the weights are super small. Just make them 1s using this command:
R = R > 0;
H = H > 0;
% transpose the R matrix 
% because their code assumes the matrix has prots as rows and goids as columns 
R = R';

addpath('/home/jeffl/git-workspace/from-lit/aptrank');
addpath('/home/jeffl/git-workspace/from-lit/aptrank/aptrank')
% Split R into Rtrain and Rtest:
rho = 0.5; % the percentage of annotations used in training.
[Rtrain,Rtest] = splitR(R,rho);

% Convert H into bi-directional:
lambda = 0.5;
dH = directH(H,lambda);

% Launch MATLAB parallel computing pool:
%matlabpool('open',12); % Set the number of cores
% NOTE: If you use MATLAB/R2013b or higher version,
% please use parpool() instead of matlabpool().
% my local machines have a limit of 4 parallel processes
parpool('local',4);
% For more information about parpool(),
% please refer to http://www.mathworks.com/help/distcomp/parpool.html

% Execution:
K = 8; % Markov chain iterations
S = 5; % Number of shuffles
t = 0.5; % Split percentage of Rtrian into Rfit and Reval
diffusion_type = 'oneway'; % Input either 'oneway' or 'twoway'.
Xa = aptrank(G,Rtrain,dH,K,S,t,diffusion_type);

% Evaluation:
auc = calcAUC(Xa,Rtrain,Rtest);
disp(['AUROC = ',num2str(auc)])
map = calcMAP(Xa,Rtrain,Rtest);
disp(['MAP = ',num2str(map)])

% Close MATLAB parallel computing pool:
%matlabpool('close');
% For MATLAB/R2013b or higher, please shut down a parallel pool using:
p = gcp; delete(p);

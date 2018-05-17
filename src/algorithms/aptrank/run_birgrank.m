function run_birgrank(network_file, annotations_and_go_dag_file)

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

% Load one of the datasets. Here take the yeast for example:
addpath('/home/jeffl/git-workspace/from-lit/aptrank');
addpath('/home/jeffl/git-workspace/from-lit/aptrank/birgrank');

% Split R into Rtrain and Rtest:
rho = 0.5; % the percentage of annotations used in training.
[Rtrain,Rtest] = splitR(R,rho);

% Convert H into bi-directional:
lambda = 0.5;
dH = directH(H,lambda);

% Execution:
%addpath('birgrank')
alpha = 0.5; % PageRank coefficient
theta = 0.5; % percentage of Rtrain in seeding
mu = 0.5; % percentage of random walks staying in G
Xh = birgrank(G,Rtrain,dH,alpha,theta,mu);

% Evaluation:
auc = calcAUC(Xh,Rtrain,Rtest);
disp(['AUROC = ',num2str(auc)])
map = calcMAP(Xh,Rtrain,Rtest);
disp(['MAP = ',num2str(map)])


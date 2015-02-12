% clear all; clc;
tic;
% Define date and run index
% Date = '2014_12_25';
% Index = 'Run1';
% Folder = strcat('./',Date,'/',Index,'/');
global Folder
%% Read files
ESig_X0 = dlmread(strcat(Folder,'ESig_X0.dat'),'\t');
ESig_X1 = dlmread(strcat(Folder,'ESig_X1.dat'),'\t');
ESig_X2 = dlmread(strcat(Folder,'ESig_X2.dat'),'\t');
ESig_X3 = dlmread(strcat(Folder,'ESig_X3.dat'),'\t');
delete(strcat(Folder,'ESig_X*.dat'));

% ESig_X = ESig_X0;
ESig_X = ESig_X0-ESig_X1+ESig_X2-ESig_X3;
toc;
dlmwrite(strcat(Folder,'ESig_X.dat'),ESig_X,'\t');
toc;
% clear all; clc;
clear all;
global Folder
% Date = '2014_12_29';
% Index = 'Run1';
% Folder = strcat('./',Date,'/',Index,'/');

Folder = './';

%Phase, wells, L, EIS, EID, gamma, E, Probe
TwoLevel_Fn_Prop(0,10,15,0,0,0.73,5E-2*pi,1)
TwoLevel_Fn_Prop(1,10,15,0,0,0.73,5E-2*pi,1)
TwoLevel_Fn_Prop(2,10,15,0,0,0.73,5E-2*pi,1)
TwoLevel_Fn_Prop(3,10,15,0,0,0.73,5E-2*pi,1)
CombineMat
PlotSpecs

% %NEXT RUN
% clear all;
% 
% Index = 'EIS1wells3';
% if ~isdir(strcat('./',Index,'/'))
%     mkdir(strcat('./',Index,'/'))
% end
% Folder = strcat('./',Index,'/');
% 
% %Phase, wells, L, EIS, EID, gamma, E
% TwoLevel_Fn_Prop(0,3,0,1,0,0.5,5E-2*pi,1)
% TwoLevel_Fn_Prop(1,3,0,1,0,0.5,5E-2*pi,1)
% TwoLevel_Fn_Prop(2,3,0,1,0,0.5,5E-2*pi,1)
% TwoLevel_Fn_Prop(3,3,0,1,0,0.5,5E-2*pi,1)
% CombineMat
% PlotSpecs
% 
% %NEXT RUN
% clear all;
% 
% Index = 'EID1gam1wells1';
% if ~isdir(strcat('./',Index,'/'))
%     mkdir(strcat('./',Index,'/'))
% end
% Folder = strcat('./',Index,'/');
% 
% %Phase, wells, L, EIS, EID, gamma, E
% TwoLevel_Fn_Prop(0,1,0,0,1,0.1,5E-2*pi,1)
% TwoLevel_Fn_Prop(1,1,0,0,1,0.1,5E-2*pi,1)
% TwoLevel_Fn_Prop(2,1,0,0,1,0.1,5E-2*pi,1)
% TwoLevel_Fn_Prop(3,1,0,0,1,0.1,5E-2*pi,1)
% CombineMat
% PlotSpecs
% 
% %NEXT RUN
% clear all;
% 
% Index = 'EID10gam1wells1';
% if ~isdir(strcat('./',Index,'/'))
%     mkdir(strcat('./',Index,'/'))
% end
% Folder = strcat('./',Index,'/');
% 
% %Phase, wells, L, EIS, EID, gamma, E
% TwoLevel_Fn_Prop(0,1,0,0,10,0.1,5E-2*pi,1)
% TwoLevel_Fn_Prop(1,1,0,0,10,0.1,5E-2*pi,1)
% TwoLevel_Fn_Prop(2,1,0,0,10,0.1,5E-2*pi,1)
% TwoLevel_Fn_Prop(3,1,0,0,10,0.1,5E-2*pi,1)
% CombineMat
% PlotSpecs
% 
% %NEXT RUN
% clear all;
% 
% Index = 'E10wells3';
% if ~isdir(strcat('./',Index,'/'))
%     mkdir(strcat('./',Index,'/'))
% end
% Folder = strcat('./',Index,'/');
% 
% %Phase, wells, L, EIS, EID, gamma, E
% TwoLevel_Fn_Prop(0,3,0,0,0,0.5,10E-2*pi,1)
% TwoLevel_Fn_Prop(1,3,0,0,0,0.5,10E-2*pi,1)
% TwoLevel_Fn_Prop(2,3,0,0,0,0.5,10E-2*pi,1)
% TwoLevel_Fn_Prop(3,3,0,0,0,0.5,10E-2*pi,1)
% CombineMat
% PlotSpecs
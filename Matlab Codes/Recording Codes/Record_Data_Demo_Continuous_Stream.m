%%==================================================================================%%
%          Streaming Data from one X4 radar - for UK DRI Launch on 18/05/2022.       %
%%==================================================================================%%
clear, clc, close all;

%% Add function paths
addpath('../../Matlab Codes/Recording Codes/Supporting Functions/');
addModuleConnectorPath();

%% Specify device Com-Port
% To find com-port easily, open XeThru Explorer
% Note: Close XeThru Explorer before running code.
COMPORT_1 = 'COM3'; % Radar COM-port

%% Call appropriate function - Done automatically.
select_recording_fun_streaming(COMPORT_1);
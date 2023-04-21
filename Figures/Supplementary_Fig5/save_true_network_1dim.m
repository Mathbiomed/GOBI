clc;
clear;
close all;


%% Fr and KF
true_network = [
    [0, -1,0];
    [0, 0,1];
    [0,0,0]];

filename = 'true_network_SFL';
save(filename, 'true_network')


%% Fr and KF
true_network = [
    [0, 1,0];
    [0, 0,1];
    [-1,0,0]];

filename = 'true_network_Fr';
save(filename, 'true_network')
filename = 'true_network_KF';
save(filename, 'true_network')


%% GW
true_network = [
    [0, 1,0,0];
    [0, 0,1,0];
    [0, 0,0,1];
    [-1,0,0,0]];

filename = 'true_network_GW';
save(filename, 'true_network')

%% Gb
true_network = [
    [0, 0,0,0,0];
    [0, 0,0,0,0];
    [0, 0,0,0,0];
    [0, 0,0,0,1];
    [-1,0,0,0,0]];

filename = 'true_network_Gb';
save(filename, 'true_network')

%% cAMP
true_network = [
    [0,0,0,0,0,1,0];
    [0,0,0,0,0,0,0];
    [0,0,0,-1,0,0,0];
    [0,0,0,0,0,0,0];
    [0,1,0,0,0,0,0];
    [0,0,0,0,0,0,1];
    [0,0,0,0,0,0,0]];

filename = 'true_network_cAMP';
save(filename, 'true_network')

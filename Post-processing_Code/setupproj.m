

% test_resultcomparison.m
% test_visualization.m
% test_plotstageheartfailure.m
% test_dnnmodelevaluation.m
output_path = '/media/maousi/Raw/lvad/dnns/dnn_6_layers_64_neurons_test/outputs/';
test_file_exact = [output_path, 'Ursino1998Model_VAD2_output_50_exact.mat'];

% For sensitivity analysis
%output_path = '/media/maousi/Raw/sensitivity/epsilon_0.2/P3_Param_LeftVentricle_kE/outputs/';
%test_file_exact = [output_path, 'Ursino1998Model_VAD2_output_88.mat'];



% Access scripts from Pre-processing/
addpath('../Pre-processing/');


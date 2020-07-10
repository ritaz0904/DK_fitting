%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Generates an option file with the specified DKI parameters
%
%   numMatPool: number of Matlab Pools to open (3)
%   dim: dimensions of file (w,h,slice,grad)
%   fname: cell containing file names (.rec uint16)
%	gradientFile: gradient file name (.txt)
%   b_value: b-value vector
%   maskname: file name of mask for skull stripping (.mat)
%   S0thresh: intensity threshold for zero-gradient image
%   gradZeroIdx: specifies which gradient index is the zero gradient
%   outputFile: name of the output file (.mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.numMatPool = 3;
options.dim = [256 256 37 32];
options.vox = [0.6 0.6 2];
options.FWHM = [2 2 2];
options.fname{1} = 'J:\MonkeyPhantom\Frenzy_080914\Frenzy_080914_9_1.REC';
options.fname{2} = 'J:\MonkeyPhantom\Frenzy_080914\Frenzy_080914_10_1.REC';
% options.fname{3} = 'recFiles/42L_b3500_4_zs_air.REC';
% options.fname{3} = 'recFiles/42L_b3500_3b_avgb0.rec';
options.gradientFile = 'E:\Tianjia\Data\Data_for_Modelling\MonkeyPhantom\Frenzy_080914\DKI\Jones30_gradient.txt';
% options.b_value = [1500 2500 3500];
options.b_value = [1500 4500];
options.maskname = '';
options.S0thresh = 50;
options.gradZeroIdx = 31;
options.outputFile = 'E:\OneDrive\Frenzy_DKIoutput.mat';
save E:\OneDrive\Frenzy_DKIoutput.mat options

%% dallas data
options.numMatPool = 6;
options.dim = [256 256 65 32];
options.vox = [224/256 224/256 2.2];
options.FWHM = [2 2 2];
options.fname{1} = 'E:\Tianjia\DKI_Tianjia\DKI_examples\3t7238\raw_data\3t7238_7_1.rec';
options.fname{2} = 'E:\Tianjia\DKI_Tianjia\DKI_examples\3t7238\raw_data\3t7238_8_1.rec';
% options.fname{3} = 'recFiles/42L_b3500_4_zs_air.REC';
% options.fname{3} = 'recFiles/42L_b3500_3b_avgb0.rec';
options.gradientFile = 'C:\Users\HuangLab\OneDrive\_DKI_Austin_fromMinhui_2018\DKIcode_v2_Austin\gradFiles\Jones30_gradient_1.txt';
options.b_value = [1000 2500];
options.maskname = '';
options.S0thresh = 10;
options.gradZeroIdx = 1;
options.outputFile = 'E:\Tianjia\DKI_Tianjia\DKI_examples\3t7238\3t7238_061819_1';
save E:\Tianjia\DKI_Tianjia\DKI_examples\3t7238\Options_3t7238_061819_1 options
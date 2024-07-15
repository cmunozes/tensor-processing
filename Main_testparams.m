clear 
close all

load('my_colourmaps.mat')
cmap_fa = [0 0 0; flipud(othercolor('RdYlGn10'))];

avg   = [0 1]; % True
sigflt = [1 3]; % 1
szflt  = [5 11 15]; % 11

% param = combvec(sigflt,szflt);

FL50um  = tiffreadVolume('/Users/cmunozes/Documents/GitHub/st_fibre_analysis_hbp/A8442_B3_T1W50um/T1W50um_stack.tif');
FL50um  = double(FL50um)/(2^16-1);

TSE65um  = tiffreadVolume('/Users/cmunozes/Documents/GitHub/st_fibre_analysis_hbp/A8442_B3_TSE65um/TSE65um_stack.tif');
TSE65um  = double(TSE65um)/(2^16-1);

TSE200um  = tiffreadVolume('/Users/cmunozes/Documents/GitHub/st_fibre_analysis_hbp/A8447_TSE200um/TSE0p2_stack.tif');
TSE200um  = double(TSE200um)/(2^8-1);
TSE200um  = TSE200um(351:1150,:,:);

inputImg = TSE200um;
%%
level = graythresh(inputImg(:,:));
mask  = inputImg >level;
mask  = imclose(mask,strel("sphere",3));
mask  = imfill(mask,"holes");

[t11, t12, t13, t22, t23, t33] = compute_structure_tensor3d(inputImg,...
    'average',true,'sizeAveragingFilter',7,'sigmaAveragingFilter',1);

[lambda1, lambda2, lambda3] = compute_eigenvalues_of_tensor3d(t11, t12, t13, t22, t23, t33);

[eigenvalues, eigenvectors] = compute_eig_of_tensor3d(t11, t12, t13, t22, t23, t33, mask);
%%
HA = zeros(size(mask));

for ss=1:size(inputImg,3)
    
    [HA(:,:,ss),long_axis_info_temp]...
        = helix_angle_approx(mask(:,:,ss), ...
        squeeze(eigenvectors(:,:,ss,:,:)),[],"LV_SA");
    
    long_axis_info{ss}.vector=long_axis_info_temp.vector;
    long_axis_info{ss}.points=long_axis_info_temp.points;
    long_axis_info{ss}.curve=long_axis_info_temp.curve;
    long_axis_info{ss}.divider=long_axis_info_temp.divider;
    
    % for some noisy pixels
    HA(:,:,ss) = HA(:,:,ss);
    
end


%% Test figures

for nn = 1:size(inputImg,3)
    figure(1),
    subplot(1,2,1), imshow(inputImg(:,:,nn).*mask(:,:,nn),[0 1],'InitialMagnification',1500)
    subplot(1,2,2), imshow(HA(:,:,nn),[-90 90],'InitialMagnification',1500), colormap(my_bl_rd_wh_bl_bl_cmap)    
    pause(0.2)
end
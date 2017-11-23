function [obj_stats,seg_stats,features] = quantify_single_section(atlas,seg,slice,...
    output_dir,obj_lbl,curr_atlas,metadata)
%QUANTIFY_SINGLE_SECTION quantify objects in a section
% Takes resolution from the original text file
% Calculates the downscale ratio using original and downscale height and
% weight
% Read the Atlas cut (ReadSlice.m), resize it to the right size
% Read the segmentation and create a binary image using the information
% from the label of the object of interest
% Remove unique point using a morphological operation (opening) with a
% diamond kernel of size 1
% Count individual connected objects
% Populate teh JSON file with properties of each individual object
% Plot the output


%%% Read the slice
slice_im = imread(slice);
slice_im_size = size(slice_im);

%%% Read the atlas segmentation and color the slice
atlas_im = ReadSlice(atlas);
atlas_orig_size = size(atlas_im);

%%% Resize the atlas to the size of the section
atlas_im = imresize(atlas_im,slice_im_size(1:2),'method','nearest');
atlas_im_size = size(atlas_im);

%%% Get the pixel area in the downsampled image (def: no info)
pixel_dim        = curr_atlas.pixel_resolution*(atlas_orig_size(2)./slice_im_size(2));
pixel_area       = pixel_dim*pixel_dim;
pixel_area_units = curr_atlas.pixel_resolution_unit;

%%% Color the bin image
[atlas_im_rgb,lbl_lst,lbl_idx,lbl_pixel,lbl_clr] = colorAtlasImages(atlas_im,curr_atlas.atlas_lbl_file);

%%% Make it a structure
for iR = 1:length(lbl_lst)
    seg_stats(iR).name  = lbl_lst{iR};
    seg_stats(iR).idx   = lbl_idx(iR);
    seg_stats(iR).rgb   = lbl_clr(iR,:);
    seg_stats(iR).pixel = lbl_pixel(iR);
    seg_stats(iR).area  = lbl_pixel(iR) * pixel_area;
    seg_stats(iR).area_units = [pixel_area_units 'x' pixel_area_units];
end

%% Read the segmentation and create a object binary file
[~,sl_name,~] = fileparts(seg);
[seg_im,seg_cmap] = imread(seg);
if size(seg_im,3)>1
    seg_im = rgb2ind(seg_im,seg_cmap,'nodither');
end
seg_im_size = size(seg_im);
obj_im = zeros(size(seg_im),'like',seg_im);
obj_im(seg_im==obj_lbl) = 1;

%% Remove small artefacts
% sent on 26/04 to collaborators
k_def = strel('diamond',1);
obj_im_cl = imopen(obj_im,k_def);
% figure('Position',[100 100 900 600]);subplot(1,2,1);imagesc(obj_im);axis image;axis off;subplot(1,2,2);imagesc(obj_im_cl);axis image;axis off

%% Check dimensions
assert(all(atlas_im_size(1:2)==slice_im_size(1:2)),...
    'Size atlas (%d %d) is different from size image (%d %d)',...
    atlas_im_size(1),atlas_im_size(2),slice_im_size(1),slice_im_size(2));
assert(all(atlas_im_size(1:2)==seg_im_size(1:2)),...
    'Size atlas (%d %d) is different from size seg (%d %d)',...
    atlas_im_size(1),atlas_im_size(2),seg_im_size(1),seg_im_size(2));

%% Returns measurements for the set of properties specified by
% properties for each connected component (object) in the binary image
stats  = regionprops(logical(obj_im_cl),'Centroid','Area','BoundingBox',...
    'Orientation', 'MajorAxisLength','MinorAxisLength');
statsR = regionprops(logical(obj_im_cl),slice_im(:,:,1),'MeanIntensity');
statsG = regionprops(logical(obj_im_cl),slice_im(:,:,2),'MeanIntensity');
statsB = regionprops(logical(obj_im_cl),slice_im(:,:,3),'MeanIntensity');
% A = [stats.Area];
n_obj = length(stats)-1;
iR = 0;
%
for iL = 1:n_obj
    %% Count the cells
    % only include it if there is an atlas label associated to the location
    % of the centroid of the object
    if atlas_im(round(stats(iL).Centroid(2)),round(stats(iL).Centroid(1)))>0
        iR = iR + 1;
        %%% Label connected regions
        % Area properties
        obj_stats(iR).object_pixel    = stats(iL).Area; %#ok<*AGROW>
        obj_stats(iR).object_area     = stats(iL).Area * pixel_area;
        obj_stats(iR).object_area_units = [pixel_area_units 'x' pixel_area_units];
        
        %Location properties (x,y)
        obj_stats(iR).object_centroid_pixel = stats(iL).Centroid;
        
        %Location in ABA space: careful with centroid coordinates (x,y) and
        %image size (height(vertical) width(horizontal)) and atlas standards
        if ~(isempty(metadata.o_vec) || isempty(metadata.u_vec) || isempty(metadata.v_vec))
        objcentpix_height_width_norm =...
            [obj_stats(iR).object_centroid_pixel(2),...
            obj_stats(iR).object_centroid_pixel(1)]./seg_im_size;
        if strcmpi(curr_atlas.orientation(1),'l')
        obj_stats(iR).object_centroid_atlas = metadata.o_vec +...
                    metadata.u_vec * (1-objcentpix_height_width_norm(2)) +...
                    metadata.v_vec * objcentpix_height_width_norm(1);
        
        else
                    obj_stats(iR).object_centroid_atlas = metadata.o_vec +...
                    metadata.u_vec * objcentpix_height_width_norm(2) +...
                    metadata.v_vec * objcentpix_height_width_norm(1);
        end
        obj_stats(iR).object_atlas_anchored = curr_atlas.name;
        end
        
        % Verification that the coordinates calculated are within the ABA
        % space size
        if ~all(obj_stats(iR).object_centroid_atlas < curr_atlas.dim)
            % Bad news, just put NaNs
            obj_stats(iR).object_centroid_atlas = [NaN NaN NaN];
        end
        
        %Shape properties
        obj_stats(iR).object_ori            = stats(iL).Orientation;
        obj_stats(iR).object_major_al_pixel = stats(iL).MajorAxisLength;
        obj_stats(iR).object_minor_al_pixel = stats(iL).MinorAxisLength;
        % Intensity properties
        obj_stats(iR).object_RGB_mean = [statsR(iL).MeanIntensity statsG(iL).MeanIntensity statsB(iL).MeanIntensity];
        % Region belonging properties
        obj_stats(iR).region_idx      = atlas_im(round(stats(iL).Centroid(2)),round(stats(iL).Centroid(1)));
        obj_stats(iR).region_name     = lbl_lst{lbl_idx==obj_stats(iR).region_idx};
        obj_stats(iR).region_rgb      = squeeze(atlas_im_rgb(round(stats(iL).Centroid(2)),round(stats(iL).Centroid(1)),:))';
        % Slice information
        obj_stats(iR).slice_name      = sl_name;
        obj_stats(iR).slice_size      = seg_im_size;
        %
        %%%%%%%%%%%%
        %create the spatial coord json file
        features(iR).type = 'Feature';
        features(iR).geometry.type = 'Point';
        features(iR).geometry.coordinates = obj_stats(iR).object_centroid_atlas;
        features(iR).geometry.space = 'WHS_SD_rat_v1.01';
        features(iR).properties.id = obj_stats(iR).slice_name;
        features(iR).properties.region_name_str = obj_stats(iR).region_name;
        %    
    end
end
% Write images in a specific folder
output_dir_img = fullfile(output_dir,'qc_fig');
if ~exist(output_dir_img,'dir')
    mkdir(output_dir_img);
end
%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the result for this slice using ellipses
phi = linspace(0,2*pi,50);
cosphi = cos(phi);
sinphi = sin(phi);
%
warning('off', 'Images:initSize:adjustingMag');
hF=figure;
imshow(slice_im,'Border','tight');
hold on
centroids = cat(1, obj_stats.object_centroid_pixel);
%
for iC=1:length(centroids)
    plot(centroids(iC,1),centroids(iC,2),'Color',obj_stats(iC).region_rgb,'Marker','x','MarkerSize',3);
    %
    xbar = centroids(iC,1);
    ybar = centroids(iC,2);
    
    a = obj_stats(iC).object_major_al_pixel/2;
    b = obj_stats(iC).object_minor_al_pixel/2;
    
    theta = pi*obj_stats(iC).object_ori/180;
    R = [ cos(theta)   sin(theta)
        -sin(theta)   cos(theta)];
    
    xy = [a*cosphi; b*sinphi];
    xy = R*xy;
    
    x = xy(1,:) + xbar;
    y = xy(2,:) + ybar;
    %
    plot(x,y,'Color',obj_stats(iC).region_rgb,'LineWidth',1);
    %
end
hold off
%
F = getframe(hF);
imwrite(F.cdata,fullfile(output_dir_img,[sl_name '_classification.png']));
close(hF);
%%%%%%%%%%%%%%%%%%%%%%%%%
% Change atlas background to white
atlas_im_r  = uint8(ones(size(atlas_im_rgb,1),size(atlas_im_rgb,2))*255);
atlas_im_g  = uint8(ones(size(atlas_im_rgb,1),size(atlas_im_rgb,2))*255);
atlas_im_b  = uint8(ones(size(atlas_im_rgb,1),size(atlas_im_rgb,2))*255);
%
idx = (atlas_im_rgb(:,:,1)==0)&(atlas_im_rgb(:,:,2)==0)&(atlas_im_rgb(:,:,3)==0);
%
tmpr = atlas_im_rgb(:,:,1);
atlas_im_r(~idx) = tmpr(~idx);
tmpg = atlas_im_rgb(:,:,2);
atlas_im_g(~idx) = tmpg(~idx);
tmpb = atlas_im_rgb(:,:,3);
atlas_im_b(~idx) = tmpb(~idx);
%
atlas_im_rgb = cat(3,atlas_im_r,atlas_im_g,atlas_im_b);
%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the atlas
hF2=figure;
imshow(atlas_im_rgb,'Border','tight');
F2 = getframe(hF2);
imwrite(F2.cdata,fullfile(output_dir_img,[sl_name '_atlasrgb.png']));
close(hF2);
%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the blend
hF3=figure;
C=imlincomb(0.4,atlas_im_rgb,0.6,slice_im);
imshow(C,'Border','tight');
F3 = getframe(hF3);
imwrite(F3.cdata,fullfile(output_dir_img,[sl_name '_blend.png']));
close(hF3);
%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot all together
%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the blend
hF4=figure;
C=imlincomb(0.4,atlas_im_rgb,0.6,slice_im);
imshow(C,'Border','tight');
hold on
centroids = cat(1, obj_stats.object_centroid_pixel);
%
for iC=1:length(centroids)
    plot(centroids(iC,1),centroids(iC,2),'Color',[0 0 0],'Marker','x','MarkerSize',3);
end
hold off
F4 = getframe(hF4);
imwrite(F4.cdata,fullfile(output_dir_img,[sl_name '_blend_objects.png']));
close(hF4);

% Prepare the output
obj_stats = obj_stats';
seg_stats = seg_stats';
features  = features';
%
return

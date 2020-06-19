close all
clear

set(0,'DefaultFigureWindowStyle','docked');
addpath(genpath('functions'));

base_path = 'C:\Users\danie\Dropbox\CCM\Projects\Electron_Microscopy\input-images';
folders = dir(base_path);
ResultTable = table();


plugin_name = 'Hydronepherosis';
thresh_low = 85;
thresh_high = 150;
debug_level = 'Result Only';
debug_level = 'All';
threshold_smooth_param = 7;
watershed_smooth_param = 15;
seeds = false;
boarder_clear = false;
boarder_clear = 0;
intensity_filter_param = 19;
% intensity_filter_param = 30;
max_area = Inf;
min_area = 100;

save_figs = false;
save_mag = '-native';

date_str = datestr(now,'yyyymmddTHHMMSS');
% save_path_prefix = sprintf('plots/%s_%s/',date_str, strrep(folder_name,'\','__'));
save_path_prefix = sprintf('%s/',date_str);
mkdir(['plots/' save_path_prefix]);


for folder = folders'
  if ~folder.isdir | strcmp(folder.name,'.') | strcmp(folder.name,'..')
    continue
  end

  folder_name = folder.name;
  full_path = fullfile(folder.folder, folder.name);
  files = dir(sprintf('%s/*.tif', full_path));
  if isempty(files)
    % instead of finding images, we found a folder in a folder, now lets look inside that folder
    files = dir(full_path);
    for file = files'  % should be '.', '..', and a directory containing images
      if file.isdir
        folder = file;
      end
    end
  end
  full_path = fullfile(folder.folder, folder.name);
  files = dir(sprintf('%s/*.tif', full_path));

  for file=files'
    close all;
    file_path = fullfile(file.folder, file.name)
    % file_path = 'C:\Users\danie\Dropbox\CCM\Projects\Electron_Microscopy\From-Claire2\im.tif'
    orig_img = imread(file_path);
    img = orig_img;
    img_size_y = size(img,1);
    img_size_x = size(img,2);

    % Skip if three dimensions because that means a weird overlay is happening
    if ndims(orig_img)==3
      warning(sprintf('Skipping image because it has a color overlay: %s', file.name));
      continue;
      % img=rgb2gray(img);
    end

    % figure;imshow(img);

    % Smooth
    % img_smooth = imgaussfilt(img,threshold_smooth_param);
    % if ismember(debug_level,{'All'})
    %   f = figure(886); clf; set(f,'name','smooth for threshold','NumberTitle', 'off');
    %   imshow(img_smooth,[]);
    % end

    %% MAIN
    % figure; imshow(img,[]);

    % Threshold by intensity
    img_thresh = img < 184 & img > 77;
    % figure; imshow(img_thresh,[]);

    % Threshold by smoothness
    img_smooth = imgaussfilt(double(img),5);
    img_grad=imgradient(img_smooth);
    img_grad_thresh = img_grad<prctile(img_grad(:),50);

    % Improve shape
    img_mask = img_thresh & img_grad_thresh;
    img_mask = bwareafilt(img_mask,[50000 Inf]);
    % figure; imshow(img_mask,[]);
    img_mask =imopen(img_mask,strel('disk',4));
    % figure; imshow(img_mask,[]);
    img_mask =imclose(img_mask,strel('disk',7));
    % figure; imshow(img_mask,[]);
    img_mask = bwareafilt(img_mask,[50000 Inf]);
    % figure; imshow(img_mask,[]);
    img_mask = imfill(img_mask,'holes');
    img_smooth = imgaussfilt(double(img_mask),15);
    img_mask = img_smooth>.2;
    % figure; imshow(img_mask,[]);
    % Remove unshapely objects. Keep only long slender objects.
    img_labelled = bwlabel(img_mask);
    stats = regionprops(img_labelled,'solidity');
    solidity = cat(1,stats.Solidity); % double check this
    solidity_threshold = 0.6;
    img_labelled(ismember(img_labelled,find(solidity > solidity_threshold)))=0;
    % figure; imshow(img_labelled,[]);
    img_mask_basement = img_labelled>0;

    %% Watershed foot processes
    img_smooth = imgaussfilt(double(img),15);
    img_hmin = imhmin(img_smooth,6);
    [seeds]=imregionalmin(img_hmin);
    img_smooth = imgaussfilt(double(img),7);
    img_min = imimposemin(img_smooth,seeds);
    img_ws = watershed(img_min);
    % figure; imshow(img,[]);
    % hold on
    % [xm,ym]=find(seeds);
    % hold on
    % plot(ym,xm,'or','markersize',2,'markerfacecolor','r','markeredgecolor','r')

    % break

    %% Threshold foot processes
    img_mask = img < 105;
    % figure; imshow(img_mask,[]);
    img_mask = bwareafilt(img_mask,[300 Inf]);
    % figure; imshow(img_mask,[]);
    img_mask =imopen(img_mask,strel('disk',4));
    % figure; imshow(img_mask,[]);
    % img_mask =imclose(img_mask,strel('disk',7));
    % figure; imshow(img_mask,[]);
    img_mask = bwareafilt(img_mask,[300 Inf]);
    % figure; imshow(img_mask,[]);
    img_smooth = imgaussfilt(double(img_mask),4);
    img_mask = img_smooth>.5;
    % img_mask = imfill(img_mask,'holes');
    % figure; imshow(img_mask,[]);
    % img_mask(img>85)=0;
    % figure; imshow(img_mask,[]);
    img_mask_foot = img_ws & img_mask;
    % figure; imshow(img_mask_foot,[]);
    img_mask_basement_eroded = imdilate(img_mask_basement,strel('disk',5));
    img_overlap = img_mask_foot & img_mask_basement_eroded;
    img_good_foot = imreconstruct(img_overlap,img_mask_foot);
    img_basement_skel = bwskel(img_mask_basement);
    img_good_foot(img_basement_skel) = 0;
    img_good_foot(img_mask_basement) = 0;
    img_good_foot = bwareafilt(img_good_foot,[600 50000]);
    figure; imshow(img_good_foot,[]);
    % img_foot_labelled = bwlabel(img_good_foot);
    % stats = regionprops('table', img_foot_labelled, img, 'MeanIntensity');
    % img_foot_labelled(ismember(img_foot_labelled,find(stats.MeanIntensity > 80)))=0;
    % img_mask_foot = img_foot_labelled>0;
    % figure; imshow(img_foot_labelled,[]);
    % figure; imshow(img,[]);

    %% Remove image perimeter from basement
    img_perim_basement = bwperim(img_mask_basement);
    img_perim_basement(1,1:end)=0;
    img_perim_basement(end,1:end)=0;
    img_perim_basement(1:end,1)=0;
    img_perim_basement(1:end,end)=0;
    % figure; imshow(img_perim_basement,[]);
    labelled_mask_basement = bwlabel(img_mask_basement);

    % Count feet for each side of each basement
    all_foot_count = [];
    all_foot_mean_int = [];
    all_basement_edge_stdev = [];
    for base_num=1:max(labelled_mask_basement(:))
      one_perim_basement = labelled_mask_basement==base_num;
      one_basement_perim = bwlabel(one_perim_basement & img_perim_basement);
      for side_num=1:max(one_basement_perim(:))
        one_basement_side = one_basement_perim==side_num;
        % figure; imshow(one_basement_side,[]);
        one_basement_side_eroded = imdilate(one_basement_side,strel('disk',5));
        img_overlap = one_basement_side_eroded & img_good_foot;
        one_side_foot = imreconstruct(img_overlap,img_good_foot);
        % figure; imshow(one_side_foot,[]);
        % Calc foot count
        foot_count = max(max(bwlabel(one_side_foot)));
        all_foot_count(base_num, side_num) = foot_count;
        % Calc foot mean pixel intensity (to resolve ties)
        foot_pixels = img(one_side_foot);
        all_foot_mean_int(base_num, side_num) = mean(foot_pixels(:));
        % Calc basement edge average variation (for eliminating uninteresting basements)
        one_basement_side_eroded = imdilate(one_basement_side,strel('disk',50));
        one_basement_side_eroded(img_mask_basement)=0;
        % figure; imshow(one_basement_side_eroded,[]);
        pixels = double(img(one_basement_side_eroded));
        all_basement_edge_stdev(base_num, side_num) = std(pixels(:));
      end
    end
    % For each basement, delete all feet on the side of the basement with fewer feet
    for base_num=1:max(labelled_mask_basement(:))
      one_perim_basement = labelled_mask_basement==base_num;
      one_basement_perim = bwlabel(one_perim_basement & img_perim_basement);
      for side_num=1:max(one_basement_perim(:))
        % Resolve ties, ie. same number of feet on both sides, remove the mean brighter side, brighter is less likely to be feet
        if size(unique(all_foot_count(base_num, :)),2)==1 && all_foot_mean_int(base_num, side_num) == max(all_foot_mean_int(base_num, :)) % Exact same number of feet on both sides
          one_basement_side = one_basement_perim==side_num;
          % figure; imshow(one_basement_side,[]);
          one_basement_side_eroded = imdilate(one_basement_side,strel('disk',5));
          img_overlap = one_basement_side_eroded & img_good_foot;
          one_side_foot = imreconstruct(img_overlap,img_good_foot);
          % do delete
          img_good_foot(one_side_foot) = 0;
          fprintf('Deleting by int\n');
        % Delete side with fewer feet
        elseif all_foot_count(base_num, side_num) ~= max(all_foot_count(base_num, :))
          one_basement_side = one_basement_perim==side_num;
          one_basement_side_eroded = imdilate(one_basement_side,strel('disk',5));
          img_overlap = one_basement_side_eroded & img_good_foot;
          one_side_foot = imreconstruct(img_overlap,img_good_foot);
          % figure; imshow(one_basement_side_eroded,[]);
          % figure; imshow(one_side_foot,[]);
          % do delete
          img_good_foot(one_side_foot) = 0;
          fprintf('Deleting by count\n');
        % Delete feet and basement if area near basement is not variable in intensity (likely to not be foot processes which are highly variable in intensity)
        elseif all_basement_edge_stdev(base_num, side_num) < 33
          one_basement_side = one_basement_perim==side_num;
          one_basement_side_eroded = imdilate(one_basement_side,strel('disk',5));
          img_overlap = one_basement_side_eroded & img_good_foot;
          one_side_foot = imreconstruct(img_overlap,img_good_foot);
          % do delete of feet
          img_good_foot(one_side_foot) = 0;
          % do delete of basement
          one_basement = imreconstruct(logical(one_basement_perim),img_mask_basement);
          figure; imshow(one_basement,[]);
          img_mask_basement(one_basement)=0;
        end
      end
    end
    figure; imshow(img_good_foot,[]);

    img_good_foot = imfill(img_good_foot,'holes');
    % img_good_foot = imclearborder(img_good_foot);

    img_labelled = bwlabel(img_good_foot);
    stats = regionprops(img_labelled,'solidity');
    solidity = cat(1,stats.Solidity); % double check this
    solidity_threshold = 0.6;
    img_labelled(ismember(img_labelled,find(solidity < solidity_threshold)))=0;
    % figure; imshow(img_labelled,[]);
    img_good_foot = img_labelled>0;


    % Visualization
    if ismember(debug_level,{'All','Result Only','Result With Seeds'})
      f = figure(743); clf; set(f,'name',[plugin_name ' Result'],'NumberTitle', 'off')
      % Display original image
      % Cast img as double, had issues with 32bit
      img8 =  img; %im2uint8(double(img));
      if min(img8(:)) < prctile(img8(:),99.5)
          min_max = [min(img8(:)) prctile(img8(:),99.5)];
      else
          min_max = [];
      end
      imshow(img8,[min_max]);
      hold on
      % % Display color overlay (grey all objects)
      % labelled_perim = imdilate(bwperim(labelled_img_all),strel('disk',1));
      % labelled_rgb = label2rgb(uint32(labelled_perim), [.7 .7 .7], [1 1 1]);
      % himage = imshow(im2uint8(labelled_rgb),[min_max]);
      % himage.AlphaData = labelled_perim*1;
      % Display color overlay (colored possible hydronephrosis)

      % labelled_perim = imdilate(bwlabel(bwperim(img_mask_basement)),strel('disk',1));
      labelled_rgb = label2rgb(uint32(img_mask_basement),[.3 1 0], [1 1 1]);
      himage = imshow(im2uint8(labelled_rgb),[min_max]);
      himage.AlphaData = img_mask_basement*.15;
      
      labelled_rgb = label2rgb(uint32(bwperim(img_mask_basement)),[.3 1 0], [1 1 1]);
      himage = imshow(im2uint8(labelled_rgb),[min_max]);
      himage.AlphaData = bwperim(img_mask_basement)*1;
      
      labelled_perim = imdilate(bwlabel(bwperim(img_good_foot)),strel('disk',1));
      labelled_rgb = label2rgb(uint32(labelled_perim), 'jet', [1 1 1], 'shuffle');
      himage = imshow(im2uint8(labelled_rgb),[min_max]);
      himage.AlphaData = labelled_perim*1;
      % if ismember(debug_level,{'All','Result With Seeds'})
      %   if ~isequal(seeds,false)
      %     % Display red dots for seeds
      %     [xm,ym]=find(seeds);
      %     hold on
      %     plot(ym,xm,'or','markersize',2,'markerfacecolor','r','markeredgecolor','r')
      %   end
      % end
      hold off
      if save_figs
        pause(0.1)
        fig_name = sprintf('plots/%s/%s_%s.png',save_path_prefix,folder_name,file.name);
        export_fig(fig_name,save_mag);
      end
    end
  end
  % break
end
% Save to disk

% msgbox('Done!')

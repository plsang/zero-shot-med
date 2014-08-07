function resize_images(org_pathname, thumb_folder, scale)
img_files = dir(org_pathname);
org_img_folder = fileparts(org_pathname);
img_files = {img_files(:).name};
num_imgs = length(img_files);
for i = 1:num_imgs
    tic;
    img = imread(fullfile(org_img_folder,img_files{i}));
    img = imresize(img,scale);
    imwrite(img,fullfile(thumb_folder,img_files{i}),'jpg');
    fprintf('\r%d/%d %.0fs', i, num_imgs, toc);
end

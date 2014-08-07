function [score,good] = ransac_verify( f1, f2, matches,fg_idx,fg_scoring,fg_wei, im1_name, im2_name, show)
% SIFT_MOSAIC Demonstrates matching two images using SIFT and RANSAC
%
%   SIFT_MOSAIC demonstrates matching two images based on SIFT
%   features and RANSAC and computing their mosaic.
%
%   SIFT_MOSAIC by itself runs the algorithm on two standard test
%   images. Use SIFT_MOSAIC(IM1,IM2) to compute the mosaic of two
%   custom images IM1 and IM2.

% AUTORIGHTS

good = true;
score = 0;


if ~exist('show','var')
    show = false;
end

if ~exist('method','var')
    method = 'Czech';%Czech
end

if fg_scoring == 1 || fg_scoring == 3
    % only check fg region and discard other matches from the bg rgion
    old_matches = matches;
    matches = matches(:,fg_idx(int32(matches(1,:))));
elseif fg_scoring == 2
    % for fg_scoring method 2, only counting fg region, so fg_wei is not necessary
    fg_wei = 1;
end
numMatches = size(matches,2) ;
if numMatches<=4 || length(unique(matches(1,:)))<=4 || length(unique(matches(2,:)))<=4 ...
      %||  length(unique(f1(1:2,unique(matches(2,:)))', 'rows'))<=4
      %||  length(unique(f2(1:2,unique(matches(2,:)))', 'rows'))<=4
    good = false;
    return;
end

if fg_scoring == 3
    % only check fg region and discard other matches from the bg rgion
    matches = old_matches;
    numMatches = size(matches,2) ;
end
X1 = [f1(1:2,matches(1,:)); ones(1,numMatches)];
X2 = [f2(1:2,matches(2,:)); ones(1,numMatches)];

% --------------------------------------------------------------------
%                                         RANSAC with homography model
% --------------------------------------------------------------------

switch method
case 'vgg'
    addpath(genpath('./MatlabFns'));
    t = .001;  % Distance threshold for deciding outliers
    if show
        tic;
    end
    [H, inliers] = ransacfithomography(X1, X2, t);
    if sum(sum(abs(H))) < 0.00001
         good = false;
        return;
    end
    if show
        fprintf('vgg %.4f, #inliers %d',toc,length(inliers));
        H
    end
    if isempty(inliers) || length(unique(matches(2,inliers)))<=4 || ...
            length(unique(f2(1:2,unique(matches(2,inliers)))', 'rows'))<=4
        good = false;
        return;
    end
    score = length(inliers);
    ok = false(1,numMatches);
    ok(inliers) = true;
case 'Czech'
    addpath(genpath('./LORANSAC/src'));
    if show
        tic;
    end
     [X1, T1] = normalise2dpts(X1);
     [X2, T2] = normalise2dpts(X2);
     TC_PAIRS = double([X1;X2]);
    if fg_scoring == 3,
        % only using pts within fg region to compute H while scoring with all pts
        TC_PAIRS = TC_PAIRS(:,fg_idx(int32(matches(1,:))));
    end
    THRESHOLD = 0.001;
    LO = int32(1);
    CONF = .95;
    INL_LIMIT = 28;
    MAX_SAMPLES = 1000; 
    [H,inliers] = loransacH(TC_PAIRS, THRESHOLD, LO, CONF, INL_LIMIT, MAX_SAMPLES);
    % in case zero H
    if sum(sum(abs(H))) < 0.00001
         good = false;
        return;
    end
    if fg_scoring == 3,
        % Calculate, in both directions, the transfered points    
        Hx1    = H*X1;
        invHx2 = H\X2;
        
        % Normalise so that the homogeneous scale parameter for all coordinates
        % is 1.
        
        X1     = hnormalise(X1);
        X2     = hnormalise(X2);     
        Hx1    = hnormalise(Hx1);
        invHx2 = hnormalise(invHx2); 
        
        d2 = sum((X1-invHx2).^2)  + sum((X2-Hx1).^2);
        inliers = find(abs(d2) < THRESHOLD);    
    end
    %[H,inliers] = loransacH(TC_PAIRS, THRESHOLD, LO, CONF);
    inliers = find(inliers);
    if show
        [~,name] = fileparts(im2_name);
        fprintf('Czech %s %.4f, #inliers %d',name,toc,length(inliers));
        H
    end
    if isempty(inliers) || length(unique(matches(2,inliers)))<=4 || ...
            length(unique(f2(1:2,unique(matches(2,inliers)))', 'rows'))<=4
        good = false;
        return;
    end
    if fg_scoring == 2,
        % only scoring pts within fg region, differ from 1 in that bg points can help to compute H
        inliers = inliers(fg_idx(int32(matches(1,inliers))));
    end
    weight=ones(size(fg_idx));
    weight(fg_idx) = fg_wei;
    score = sum(weight(int32(matches(1,inliers))));
    ok = false(1,numMatches);
    ok(inliers) = true;
case 'vlfeat'
    iter = 100;
    H = cell(1,iter);
    ok = cell(1,iter);
    score = zeros(1,iter);
    %[x1,x2,y1,y2]=get_bounding_rect(f1);
    for t = 1:iter
        % estimate homograpyh
        subset = vl_colsubset(1:numMatches, 4) ;
        A = [] ;
        for i = subset
            A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
        end
        [U,S,V] = svd(A) ;
        H{t} = reshape(V(:,9),3,3) ;
        
        % score homography
        X2_ = H{t} * X1 ;
        du = X2_(1,:)./X2_(3,:) - X2(1,:);
        dv = X2_(2,:)./X2_(3,:) - X2(2,:);
        ok{t} = abs(du)<4 & abs(dv)<4 ;
        score(t) = sum(ok{t}) ;
        %     if score(t) > good_num
        %         score = score(t);
        %         return;
        %     end
    end
    [score, best] = max(score) ;
    H = H{best} ;
    ok = ok{best} ;
end

function err = residual(H)
u = H(1) * X1(1,ok) + H(4) * X1(2,ok) + H(7) ;
v = H(2) * X1(1,ok) + H(5) * X1(2,ok) + H(8) ;
d = H(3) * X1(1,ok) + H(6) * X1(2,ok) + 1 ;
du = X2(1,ok) - u ./ d ;
dv = X2(2,ok) - v ./ d ;
err = sum(du.*du + dv.*dv) ;
end
    

if show    
    % --------------------------------------------------------------------
    %                                                  Optional refinement
    % --------------------------------------------------------------------
    
    if exist('fminsearch') == 2
        H = H / H(3,3) ;
        opts = optimset('Display', 'none', 'TolFun', 1e-8, 'TolX', 1e-8) ;
        H(1:8) = fminsearch(@residual, double(H(1:8)'), opts) ;
    else
        warning('Refinement disabled as fminsearch was not found.') ;
    end
    
    % --------------------------------------------------------------------
    %                                                         Show matches
    % --------------------------------------------------------------------
    if ~exist('im1','var') || ~exist('im2','var')
        im1 = imread(im1_name) ;
        im2 = imread(im2_name) ;
        
        % make single
        im1 = im2single(im1) ;
        im2 = im2single(im2) ;
    end
    
    dh1 = max(size(im2,1)-size(im1,1),0) ;
    dh2 = max(size(im1,1)-size(im2,1),0) ;
    
    figure(1) ; clf ;
    subplot(2,1,1) ;
    imagesc([padarray(im1,dh1,'post') padarray(im2,dh2,'post')]) ;
    o = size(im1,2) ;
    line([f1(1,matches(1,:));f2(1,matches(2,:))+o], ...
        [f1(2,matches(1,:));f2(2,matches(2,:))]) ;
    title(sprintf('%d tentative matches', numMatches)) ;
    axis image off ;
    
    subplot(2,1,2) ;
    imagesc([padarray(im1,dh1,'post') padarray(im2,dh2,'post')]) ;
    o = size(im1,2) ;
    line([f1(1,matches(1,ok));f2(1,matches(2,ok))+o], ...
        [f1(2,matches(1,ok));f2(2,matches(2,ok))]) ;
    title(sprintf('%d (%.2f%%) inliner matches out of %d', ...
        sum(ok), ...
        100*sum(ok)/numMatches, ...
        numMatches)) ;
    axis image off ;
    
    drawnow ; 
    
    if false
        % --------------------------------------------------------------------
        %                                                               Mosaic
        % --------------------------------------------------------------------
        
        box2 = [1  size(im2,2) size(im2,2)  1 ;
            1  1           size(im2,1)  size(im2,1) ;
            1  1           1            1 ] ;
        box2_ = inv(H) * box2 ;
        box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
        box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
        ur = min([1 box2_(1,:)]):max([size(im1,2) box2_(1,:)]) ;
        vr = min([1 box2_(2,:)]):max([size(im1,1) box2_(2,:)]) ;
        
        [u,v] = meshgrid(ur,vr) ;
        im1_ = vl_imwbackward(im2double(im1),u,v) ;
        
        z_ = H(3,1) * u + H(3,2) * v + H(3,3) ;
        u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
        v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
        im2_ = vl_imwbackward(im2double(im2),u_,v_) ;
        
        mass = ~isnan(im1_) + ~isnan(im2_) ;
        im1_(isnan(im1_)) = 0 ;
        im2_(isnan(im2_)) = 0 ;
        mosaic = (im1_ + im2_) ./ mass ;
        
        figure(2) ; clf ;
        imagesc(mosaic) ; axis image off ;
        title('Mosaic') ;
        
        if nargout == 0, clear mosaic ; end
    end
    
end
end



function [u1,u2,v1,v2]=get_bounding_rect(kp)
u1=min(kp(1,:));
u2=max(kp(1,:));
v1=min(kp(2,:));
v2=max(kp(2,:));
end

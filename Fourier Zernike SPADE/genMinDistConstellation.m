% RANDOM ROOT CHAIN LINK IMPLEMENTATION
function xy = genMinDistConstellation(b, min_sep, centroid_aligned)
    % Generates the coordinates of a constellation of n_src point-sources 
    % located in a disk of radius .5 wherein each source has at least one
    % nearest neighbor that is min_sep away.
    % If align_centroid is true, the returned coordinates for a  
    % constellation whos centroid (center-of-mass) lies at the origin (0,0).
    % ------------------
    % INPUTS:
    % ------------------
    % b              : nx1 vector containing relative point weights
    % min_sep        : minimum separation distance in constellation
    % align_centroid : require centroid to be at coordinate origin
    % ------------------
    % OUTPUTS
    % ------------------
    % xy             : nx2 matrix containing the xy coordinate pairs of the
    %                  constellation. 
    %                       xy(:,1) is the x coordinate
    %                       xy(:,2) is the y coordinate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author(s): Nico Deshler, University of Arizona
    % Affiliation(s): Wyant College of Optical Sciences, University of Arizona
    % Date: March 7, 2024
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    n = numel(b);
    assert((abs(sum(b)- 1) < 1e-10) && all(b>0));
    assert( 1<n && n <= 20);
    
    D = 1;          % Diameter of support disk
    R = D/2;        % Radius of support disk
    

    %packing fractions for circle-within-a-circle up to 20 circles
    % https://en.wikipedia.org/wiki/Circle_packing_in_a_circle
    circ_pack_frac = [1,0.5000,0.6466,0.6864,0.6854,0.6666,0.7777,0.7328,0.6895,0.6878,0.7148,0.7392,0.7245,0.7474,0.7339,0.7512,0.7403,0.7609,0.8034,0.7623];

    % no possible samples if the area of the area from the optimal packing fraction
    % is less than the area that the sources may be.
    assert(n*(min_sep/2).^2 <= circ_pack_frac(n) *  (R+min_sep/2)^2); 
    
    % uniformly sample a point anywhere within the disk  of radius R-min_sep 
    % (ensures that the second point is guaranteed to lie within the disk
    % of radius R)
    p1 = sampleDiskUniform(1,R-min_sep); 
    
    % add p1 to the list of point coordinates
    xy(1,:) = p1;    
        
    % variablility in deviation from min_sep
    epsilon = min_sep/100;


    make_video = 0;
    if make_video

        
        figure('units','normalized','outerposition',[0 0 1 1])
        [th,r] = meshgrid(2*pi*linspace(0,1,100),ones(1,100)*min_sep);
        [cx,cy] = pol2cart(th(:),r(:));
        v = VideoWriter('Chain-Link','MPEG-4');
        v.Quality = 100;
        v.FrameRate = 2;
        open(v)
        writeFrame = @(v) writeVideo(v,getframe(gcf));


        set(gcf,'Color','k')
        set(gcf, 'InvertHardCopy', 'off'); 
        set(gcf,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color
        set(gcf,'GraphicsSmoothing','on')
        hold on
        scatter(xy(1,1),xy(1,2),'filled','white')
        ax = gca;
        set(ax,'Color','k')
        ax.XColor = 'w';
        ax.YColor = 'w';
        xlabel('x [rl]')
        ylabel('y [rl]')
        xlim([-.5,.5])
        ylim([-.5,.5])
        axis square
        writeFrame(v)
        plot(cx+p1(:,1),cy+p1(:,2),'Color','blue','LineWidth',1);
        writeFrame(v)
    end
    
    % graph network repn
    G = zeros(n);

    % generate remaining samples
    for k = 2:n
        % check if all the points are within the min separation criteria,
        % otherwise regenerate the kth point
        while size(xy,1) < k || ~all(pdist(xy(1:k,:)) >= min_sep - (epsilon/2))
            % sample a point on a fuzzy ring of width epsilon with radius min_sep
            [rkx,rky] = pol2cart(2*pi*rand(1), epsilon*(rand(1)-.5) + min_sep); 
            rk = [rkx,rky];
            
            % randomly pick a point to add the sampled ring point coordinates
            j = randi(k-1);
            pj = xy(j,:);   
            
            % generate the new point pk
            pk = pj + rk;
            
            % add pk to the list of point coordinates after the random
            % index j
            xy(k,:) = pk;  

           if all(pdist(xy(1:k,:)) >= min_sep - (epsilon/2))
               
                G(j,k) = 1;
           end

           if make_video
                scatter(xy(1:k,1),xy(1:k,2),'filled','white')
                %plot(graph(G(1:k,1:k),'upper'),'EdgeColor','w','NodeColor','w','NodeLabel',{},'XData',xy(:,1),'YData',xy(:,2),'Linewidth',3)
                
                if all(pdist(xy(1:k,:)) >= min_sep - (epsilon/2))
                    plot(xy([j,k],1),xy([j,k],2),'Color','w','LineWidth',1)
                else
                    plot(xy([j,k],1),xy([j,k],2),'Color','red','LineWidth',1)
                end

                writeFrame(v)
           end
        end
        if make_video
            plot(xy([j,k],1),xy([j,k],2),'Color','w','LineWidth',1)
            plot(graph(G(1:k,1:k),'upper'),'EdgeColor','w','NodeColor','w','NodeLabel',{},'XData',xy(:,1),'YData',xy(:,2),'Linewidth',3)
            plot(cx+pk(:,1),cy+pk(:,2),'Color','blue','LineWidth',1)
            writeFrame(v)
        end
    end
    
    % realign the centroid
    if centroid_aligned
        xy = xy - sum(b.*xy,1);
        
        % check if the scene still falls inside the FOV. Otherwise rerun
        % the function
        if any( sum(xy.^2,2) > R^2)
            xy = genMinDistConstellation(b, min_sep, centroid_aligned);
        end
    end
    
    if make_video
        close(v)
    end

end



function xy = sampleDiskUniform(n,R)
   % generates n samples uniformly over the unit disk
   r = R*sqrt(rand(n,1));
   th = 2*pi*rand(n,1);
   
   [x,y] = pol2cart(th,r);
   xy = [x,y];
end
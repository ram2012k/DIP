clear all
image = imread('Assign3_imgs/rectangle3.jpg');
im2 = rgb2gray(image);
edgeImage = edge(im2,'Canny');
[m,n] = size(edgeImage);
p=floor(sqrt((m-1)^2+(n-1)^2));
acc =zeros(2*p +1,2*90);
non_zeroPizel = zeros(nnz(edgeImage),2);
totalEdgePixel = nnz(edgeImage);
counter = 1
for j=1:n
    for i=1:m
        if edgeImage(i,j) == 1
            non_zeroPizel(counter,1)=i;
            non_zeroPizel(counter,2)=j;
            counter = counter+1;
        end
    end
end
counter
for i=1:totalEdgePixel
    x = non_zeroPizel(i,:);
    for j = 1:2*90
        M = (x(1)-1)*cos((j-1)*pi/180) + (x(2)-1)*sin((j-1)*pi/180);
        M =floor(M)+1;
        acc(M+p+1,j) = acc(M+p+1,j)+1;
        
    end
end


figure , imshow(image),hold on;
peaks = houghpeaks(acc,100,'Threshold',0.3*max(max(acc)));
totalP = size(peaks);
edges = zeros(m,n,'uint8');
horizontal = []
vertical = []
for i=1:totalP(1)
    
    peak = peaks(i,:);
    if peak(2) == 1
        edges(peak(1)-p-1,:) = 1;
        horizontal = [horizontal;peak];
       % plot([1,n],[peak(1)-p-1,peak(1)-p-1],'LineWidth',1,'Color','green');
    else
        if peak(2) == 91 
            edges(:,peak(1)-p-1)=1;
            vertical = [vertical;peak];
          % plot([peak(1)-p-1,peak(1)-p-1],[1,m],'LineWidth',1,'Color','green');
        else
            counter = 0;
           % acc(peak(1),peak(2));
            xmax=1 ;
            xmin=100000 ;
            ymax=0;
            ymin = 100000;
            
                diff1 = peak(2)-90;
                diff2 = peak(2)+90;
                length1 = length(find(peaks(:,2) <= diff1+2 & peaks(:,2) >= diff1-2));
                length2 = length(find(peaks(:,2) <= diff2+2 & peaks(:,2) >= diff2-2));
                if (length1 == 0 && length2 == 0)
                    disp('gdfgfd  gf ');
                    continue
                    
                end
                
                for x=1:m      
                y = (cos((peak(2)-1)*pi/180)*(1-x)/sin((peak(2)-1)*pi/180)) + ((peak(1)-p-1)/sin((peak(2)-1)*pi/180));
                y=floor(y) + 1 ;
                
               if  y > 0 && y<= n
                   if y < ymin
                       ymin=y;
                       xmin=x;
                   end
                   if y >ymax
                       ymax = y;
                       xmax=x;
                   end
                       
                       
                   counter = counter  +1 ;
                   
                   edges(x,y)=1;
               end
                end
               [xmin,xmax]
               [ymin,ymax]
            %  plot([ymin,ymax],[xmin,xmax],'LineWidth',1,'Color','green');
            
            
        end
    end
    
end
horizontal = sort(horizontal);
vertical = sort(vertical);
for i=1:length(horizontal)
    line1 = horizontal(i);
    
    line2 = vertical(1);
    pointx = line1(1)-p-1;
    point1y = line2(1)-p-1;
    
    for j=2:length(vertical)
        line2 = vertical(j);
        point2y = line2(1)-p-1;
        counter = nnz(edgeImage(pointx-1:pointx+1,point1y:point2y));
        
        if counter/(point2y-point1y) > 0.5
            plot([point1y,point2y],[pointx,pointx],'LineWidth',2,'Color','red');
        end
        point1y=point2y
    end
        
   
end

for i=1:length(vertical)
    line1 = vertical(i);
    
    line2 = horizontal(1);
    pointy = line1(1)-p-1;
    point1x = line2(1)-p-1;
    
    for j=2:length(horizontal)
        line2 = horizontal(j);
        point2x = line2(1)-p-1;
        counter = nnz(edgeImage(point1x:point2x,pointy-3:pointy+3));
        
        if counter/(point2x-point1x) > 0.5
            plot([pointy,pointy],[point1x,point2x],'LineWidth',2,'Color','red');
        end
        point1x=point2x
    end
        
   
end
            
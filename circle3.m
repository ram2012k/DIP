clear all
image = imread('Assign3_imgs/circle5.jpg');
image = imresize(image,.5);
im2 = rgb2gray(image);
edgeImage = edge(im2,'Canny');
[m,n] = size(edgeImage);
%imshow(edgeImage,[]);
rmax = floor(max(m,n)/2);
acc = zeros(m,n,rmax,'uint16');
[nrow ,ncol] = find(edgeImage);

for i=1:length(nrow)
    for r=10:rmax
        for theta=1:4:360
            a = nrow(i)-1 - r*cos(theta*pi/180);
            b = ncol(i)-1 - r*sin(theta*pi/180);
            a = floor(a) + 1;
            b = floor(b) + 1;
            if a > 0 && a <=m && b > 0 && b<=n
            
            acc(a,b,r) = acc(a,b,r)+1;
            end
        end
        
            
    end
    i
end

maximum = max(max(max(acc)));
idx = find(acc == maximum);
[v1 , v2 , v3  ] = ind2sub(size(acc),idx); 
rr = [];
cc = [];
ra = [];
neighbourhoodSize = 5;
for i=1:length(v1)
    if acc(v1(i),v2(i),v3(i)) ==0
        continue
    end
    rr = [rr;v1(i)];
    cc = [cc;v2(i)];
    ra = [ra;v3(i)];
    l1 = v1(i)-neighbourhoodSize;
    l2 = v1(i)+neighbourhoodSize;
    b1 = v2(i)-neighbourhoodSize;
    b2 = v2(i)+neighbourhoodSize;
    h1 = v3(i)-neighbourhoodSize;
    h2 = v3(i)+neighbourhoodSize;
    if l1 <=0
        l1 = 1;
    end
    if l2 > m
        l2 = m;
    end
    if b1 <= 0
        b1 = 1;
    end
    if b2 > n
        b2 = n;
    end
    if h1 <= 0
        h1 = 1;
    end
    if h2 > rmax
        h2 = rmax;
    end
    acc(l1:l2,b1:b2,h1:h2) = 0;
    
end
%non-maxima suppression
while 1
    maxi = max(max(max(acc)));
    maxi
    if maxi < 0.7*maximum
        break
    end
    idx = find(acc == maxi);
    [v1 , v2 , v3  ] = ind2sub(size(acc),idx);
    for i=1:length(v1)
        if acc(v1(i),v2(i),v3(i)) ==0
            continue
        end
        rr = [rr;v1(i)];
        cc = [cc;v2(i)];
        ra = [ra;v3(i)];
        l1 = v1(i)-neighbourhoodSize;
        l2 = v1(i)+neighbourhoodSize;
        b1 = v2(i)-neighbourhoodSize;
        b2 = v2(i)+neighbourhoodSize;
        h1 = v3(i)-neighbourhoodSize;
        h2 = v3(i)+neighbourhoodSize;
        if l1 <=0
            l1 = 1;
        end
        if l2 > m
            l2 = m;
        end
        if b1 <= 0
            b1 = 1;
        end
        if b2 > n
            b2 = n;
        end
        if h1 <= 0
            h1 = 1;
        end
        if h2 > rmax
            h2 = rmax;
        end
        acc(l1:l2,b1:b2,h1:h2) = 0;
    
    end
end
%
figure,imshow(image),hold on

%[rr , cc , ra ] = ind2sub(size(acc),find(acc > 0.66*maximum));
th = 0:pi/180:2*pi;
for i=1:length(ra)
    xunit = round((ra(i))*cos(th) + rr(i));
    yunit = round((ra(i))*sin(th) + cc(i));
    counter = 0;
    counter2=0
    for j = 1:length(xunit)
        if xunit(j) > 0 && xunit(j) <= m && yunit(j) > 0 && yunit(j) <= n 
            counter2=counter2+1;
            if edgeImage(xunit(j),yunit(j)) ==1
                counter = counter +1;
            end
        end
    end
    length(xunit)
    if counter/counter2 > 0.25
        plot(yunit,xunit,'Color','green');
    end
end

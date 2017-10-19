


Irgb = imread('Road_image_1.png');

for i = 1:size(Irgb, 1)
    for j = 1:size(Irgb, 2)
        if Irgb(i, j, 1) == 0 && Irgb(i, j, 3) == 0 
            Irgb(i,j) = 0;
        end
        if  Irgb(i, j, 2) == 0 && Irgb(i, j, 3) == 0 
            Irgb(i,j) = 0;
        end
        if  Irgb(i, j, 1) == 0 && Irgb(i, j, 2) == 0 
            Irgb(i,j) = 0;
        end
        
    end 
end



Igray = 0.2989*Irgb(:,:,1)+0.5870*Irgb(:,:,2)+0.1140*Irgb(:,:,3);

% figure
% image(Igray,'CDataMapping','scaled');
% colormap('gray')
% title('Input Image in Grayscale')
median_filt = medfilt2(Igray, [2 2]);
%imshow(median_filt);

BW = edge(Igray,'Canny');

%imshow(BW);

% for i = 1:size(BW, 1)
%     for j = 1:size(BW, 2)
%         if BW(i, j) == 0
%             BW(i, j) = 1;
%         elseif BW(i, j) == 1
%             BW(i, j) = 0;
%             
%         end
%         
%     end 
% end
copy = zeros(size(BW, 1), size(BW, 2));
num_cols = size(BW, 2); 
num_rows = size(BW, 1);
for x = 1:num_cols
    for y = 1:num_rows
        n = 0;
        if BW(y, x) == 1
            continue;
        else
            
            stop = 0;
            round = 1;
            index = 1;
            full = 1;
            while stop ~= 1
                if x-index<1 || x+index>2000 || y-index<1 || y+index>1200
                    break;
                %end
                
                elseif BW(y+index, x) == 0 && BW(y-index, x) == 0 && BW(y, x+index) == 0 && BW(y, x-index) == 0
                    n = n+1;
                    full = 0;
                else
                    stop = 1;
                    break;
                end
                %round = 2 * round - 1;
                
                pass = 0;
                for m = 1:round
                    if x-round<1 || x+round>2000 || y-round<1 || y+round>1200 
                        break;
                    end
                   % end
                    if BW(y+index, x-m) == 0 && BW(y+index, x+m) == 0 && BW(y-index, x-m) == 0 && BW(y-index, x+m) == 0 && BW(y-m, x+index) == 0 && BW(y+m, x+index) == 0 && BW(y-m, x-index) == 0 && BW(y+m, x-index) == 0
                        pass = 1;
                    else
                        pass = 0;
                        stop = 1;
                        break;
                    end
                    
                end
                if pass == 1
                        index = index + 1;
                        round = round + 1;
                        full = 1;
                end
                
            end
            if full == 1
                final = n + 0.5;
            elseif full == 0
                final = ((2*n-1)*sqrt(2))/2;
            end
            copy(y, x) = final;
        end
    end
end



 for i = 1:size(copy, 1)
     for j = 1:size(copy, 2)
         if copy(i, j) < 17
             copy(i, j) = 0;
         
         end
     end
 end
 
 
 

%  for i = 1:size(copy, 1)
%      for j = 1:size(copy, 2)
%          if copy(i, j) >= 3.5
%              copy(i, j) = 0;
%          else
%              copy(i, j) = 1;
%          end
%      end
%  end

centers = zeros(size(copy, 1), size(copy, 2));
% 
for i = 2:1:size(copy, 1)-2
     for j = 2:1:size(copy, 2)-2
         if copy(i, j) ~= 0
             max = 0;
             max_r = i;
             max_c = j;
             for r = i-1:i+1
                for c = j-1:j+1
                     if copy(r, c) > max
                         max = copy(r, c);
                         max_r = r;
                         max_c = c;
                     end
                end
             end
             if max_r==i && max_c==j
                 centers(i,j) = 1;
             end
         end
         
     end
         
end
        
% imshow(centers);        

% D1 = bwdist(BW,'euclidean');
% RGB1 = repmat(mat2gray(D1), [1 1 3]);
% imshow(RGB1);



% for i = 1:size(centers, 1)
%     for j = 1:size(centers, 2)
%         if centers(i, j) == 1
%             centers(i, j) = 0;
%         elseif centers(i, j) == 0
%             centers(i, j) = 1;
%             
%         end
%         
%     end 
% end
%imshow(max_matrix);
     
% [H,theta,rho] = hough(centers);
% figure
% imshow(imadjust(mat2gray(H)),[],...
%        'XData',theta,...
%        'YData',rho,...
%        'InitialMagnification','fit');
% xlabel('\theta (degrees)')
% ylabel('\rho')
% axis on
% axis normal 
% hold on
% colormap(gca,hot)
% % find peaks
% P = houghpeaks(H,10);
% x = theta(P(:,2));
% y = rho(P(:,1));
% plot(x,y,'s','color','black');
% lines = houghlines(BW,theta,rho,P,'FillGap',5,'MinLength',7);
% 
% figure, imshow(centers), hold on
% max_len = 0;
% for k = 1:length(lines)
%    xy = [lines(k).point1; lines(k).point2];
%    plot(xy(:,1),xy(:,2),'LineWidth',6,'Color','green');
% 
%    % Plot beginnings and ends of lines
%    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% 
%    % Determine the endpoints of the longest line segment
%    len = norm(lines(k).point1 - lines(k).point2);
%    if ( len > max_len)
%       max_len = len;
%       xy_long = xy;
%    end
% end
% % highlight the longest line segment
% plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');

% MEAN FILTER
% h = 1/3*ones(3,1);
% H = h*h';
%     % im be your image
% imfilt = filter2(H,centers);
% imshow(imfilt);


%DISTANCE TRANFORM
%D1 = bwdist(copy,'euclidean');
% %RGB1 = repmat(mat2gray(D1), [1 1 3]);
% RGB1 = repmat(mat2gray(D1));
 %imshow(D1);

% test = edge(RGB1,'Canny');
% imshow(test);
%        
 


for i = 21:1:size(centers, 1)-21
     for j = 21:1:size(centers, 2)-21
         if centers(i, j) ~= 0
             pass = 0;
             for r = i-20:i+20
                for c = j-20:j+20
                     if centers(r, c) ~= 0
                         pass = pass + 1;
                     end
                end
             end
             if pass == 1 
                 centers(i,j) = 0;
             end
         end
         
     end
         
end

imshow(centers);

s = struct('r',{},'c',{},'present',{}, 'slope', {});
for i = 1:size(centers, 1)
     for j = 1:size(centers, 2)
         index = 0;
         if centers(i, j) ~= 0
             index = index + 1;
             s(index).r = i;
             s(index).c = j;
             s(index).present = 1;
             s(index).slope;
         end
         
     end
         
end


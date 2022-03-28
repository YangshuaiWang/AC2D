 % function I = test_fcn_square(edge)
% first version use edge as input only.
% I is a vector stores upper-right vertex of square
function I = get_square(edge)
if nargin == 0;
    test_get_square();
    return;
end
% compute k to decide increasing direction, c.w. 0 then c.w. 1
edge = sortrows(edge')';
x = edge(1,:); y = edge(2,:);
k = (y(2) - y(1))/(x(2)-x(1));
x = x + eps; y = y + eps;
I = zeros(2,0);
if k == 0;
    Ix = ceil(x(1)):ceil(x(2));
    Iy = (ceil(y(1)))*ones(1,length(Ix));
%     I = sub2ind([n,m],Ix,Iy);
    I = [I, [Ix; Iy]];
    return;
end
if k == Inf;
    Iy = ceil(min(y)):ceil(max(y));
    Ix = (ceil(x(1)))*ones(1,length(Iy));
%     I = sub2ind([n,m],Ix,Iy);
    I = [I, [Ix; Iy]];
    return;
end
if k > 0
    if k > 1
        Yi = ceil(y(1))-1:ceil(y(2));
        Xi = interp1(y, x, Yi);
        Xi(1) = x(1); Xi(length(Yi)) = x(2);
        a4c = diff(ceil(Xi));
        Ix = ceil(x(1)); Iy = ceil(y(1));
%         I = zeros(1,0);
        for p = a4c
%             I = [I,sub2ind([n,m],Ix,Iy)];
            I = [I, [Ix; Iy]];
            if p == 0
                Iy = Iy + 1;
            else
                Ix = Ix + 1;
%                 I = [I,sub2ind([n,m],Ix,Iy)];
                I = [I, [Ix; Iy]];
                Iy = Iy +1;
            end
        end
    return
    else
        Xi = ceil(x(1))-1:ceil(x(2));
        Yi = interp1(x, y, Xi);
        Yi(1) = y(1); Yi(length(Xi)) = y(2);
        a4c = diff(ceil(Yi));
        Ix = ceil(x(1)); Iy = ceil(y(1));
%         I = zeros(1,0);
        for p = a4c
%             I = [I,sub2ind([n,m],Ix,Iy)];
            I = [I, [Ix; Iy]];
            if p == 0
                Ix = Ix + 1;
            else
                Iy = Iy + 1;
%                 I = [I,sub2ind([n,m],Ix,Iy)];
                I = [I, [Ix; Iy]];
                Ix = Ix +1;
            end
        end
    return
    end
end
if k < 0
    if k < -1;
        Yi = ceil(y(2))-1:ceil(y(1));
        Xi = interp1(y, x, Yi);
        Xi(1) = x(2); Xi(length(Yi)) = x(1);
        a4c = diff(ceil(Xi));
        Ix = ceil(x(2)); Iy = ceil(y(2));
%         I = zeros(1,0);
        for p = a4c
%             I = [I,sub2ind([n,m],Ix,Iy)];
            I = [I, [Ix; Iy]];
            if p == 0
                Iy = Iy + 1;
            else
                Ix = Ix - 1;
%                 I = [I,sub2ind([n,m],Ix,Iy)];
                I = [I, [Ix; Iy]];
                Iy = Iy +1;
            end
        end
        return
    else
        Xi = ceil(x(2)):-1:ceil(x(1))-1;
        Yi = interp1(x, y, Xi);
        Yi(1) = y(2); Yi(length(Xi)) = y(1);
        a4c = diff(ceil(Yi));
        Ix = ceil(x(2)); Iy = ceil(y(2));
%         I = zeros(1,0);
        for p = a4c
%             I = [I,sub2ind([n,m],Ix,Iy)];
            I = [I, [Ix; Iy]];
            if p == 0
                Ix = Ix - 1;
            else
                Iy = Iy + 1;
%                 I = [I,sub2ind([n,m],Ix,Iy)];
                I = [I, [Ix; Iy]];
                Ix = Ix - 1;
            end
        end
    end
        return
end


end
%% TEST ROUTINE
function test_get_square()
% x1 = [1.1, 7.3, 3.5]; y1 = [2, 0.8, 10.3];
% x1 = [0.2, 7.3, 3.5]; y1 = [0.3, 1.8, 10.3];
x1 = [39.4 26.205 29.3]; y1 = [-99 -88.6833 -99];
x = [39.4 26.205 29.3]; y = [-99 -88.6833 -99];
% x1 = [3, 1, 7]; y1 = [4, 6, 8];
% x1 = [3, 5, 7]; y1 = [4, 1, 8];
% x1 = [-3, -2, -2]; y1 = [4, 3, 4];
% x1 = [-5.5575, -5.5575, -4]; y1 = [3.0383, 1.5192, 1];
% x1 = [3.0383, 1.5192, 1]; y1 = [-5.5575, -5.5575, -4]; 
clf
% plot(x1([1:3,1]),y1([1:3,1]),'k');
xt = floor(min(x1)); yt = floor(min(y1));
x1 = x1-xt; y1 = y1-yt;

% plot(x1([1:3,1]),y1([1:3,1]),'k');
xmin = floor(min(x)); ymin = floor(min(y));
xmax = ceil(max(x)); ymax = ceil(max(y));
X = xmin:xmax; Y = ymin:ymax;
[XX, YY] = meshgrid(X, Y);
hold on;
plot(XX,YY,'ob')
axis([xmin-1, xmax+1, ymin-1, ymax+1]);
[m, n] = size(XX);
for i = 1:n-1
    for j = 1:m-1
        plot(X([i,i+1]),Y([j+1,j]),'--');
    end
end
plot(X,YY','b:'); plot(XX,Y,'b:');
hold off;
I3 = get_square([x1([1,2]);y1([1,2])]);
I2 = get_square([x1([1,3]);y1([1,3])]);
I1 = get_square([x1([3,2]);y1([3,2])]);
I = unique([I1, I2, I3]', 'rows')';
I = I + [xt*ones(1,length(I)); yt*ones(1,length(I))];
hold on;
% for idx = 1:length(I)
%     i = I(1,idx); j = I(2,idx);
%     xf = X([i-1,i,i,i-1]); yf = Y([j-1,j-1,j,j]);
%     fill(xf, yf, 'g');
% end
for idx = 1:length(I)
    i = I(1,idx); j = I(2,idx);
    xf = [i-1,i,i,i-1]; yf = [j-1,j-1,j,j];
    fill(xf, yf, 'g');
end
plot(x([1:3,1]),y([1:3,1]),'k');
hold off;
end

clear all;
M = dlmread('out');
n = size(M);
nf = n(1)
nc = n(2)
rows = 11;
for iter = 1:2:(nf/20)-10
  ui = M(1+rows*iter-rows:rows*iter, 1:rows);
  vi = M(rows*iter+1:rows*iter+1 + rows, 1:rows);
  idx = 1;
  for i=1:rows
      x = i;
      for j=1:nc
          y=j;
              %create the arrays for the quiver3 plot.
              xs(idx) = x;
              ys(idx) = y;
              qux(idx) = ui(i,j);
              quy(idx) = vi(i,j);
              %create matrix for the heatmap plot.
              hm(j,i) = sqrt(qux(idx)*qux(idx) + quy(idx)*quy(idx));

              idx = idx+1;
      end
  end
   iter
   pause(0.3)

%quiver(xs, ys, qux, quy, 5);
surf(hm);
title('Function plot')
xlabel('x')
ylabel('y')

%figure
hold on 
drawnow
hold off


end
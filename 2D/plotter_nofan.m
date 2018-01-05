clear all;
data = dlmread('./out/UV.txt');

n = size(data);
 data_rows = n(1);
 data_cols = n(2);
 mat_rows = 41;
 mat_cols = data_cols;
 skip = 1;
 elementSkip = 1;
 
 sys_rows = data_cols;
 sys_cols = data_cols;
 xc = sys_cols/2;
 yc = sys_rows/2; %revisar

  
pause(2);
 for base=1:skip*mat_rows*2:data_rows
     U = data(base:base+mat_rows-1, 1:mat_cols);
     V = data(base+mat_rows: base+2*mat_rows-1, 1:mat_cols);
     base
   idx = 1;
   for i=1:mat_rows
       x = i;
       for j=1:data_cols
           y=j;
               %create the arrays for the quiver3 plot.
               xs(idx) = x;
               ys(idx) = y;
               qux(idx) = U(i,j);
               quy(idx) = V(i,j);
               %create matrix for the heatmap plot.
               hm(j,i) = sqrt(qux(idx)*qux(idx) + quy(idx)*quy(idx));
 
               idx = idx+1;
       end
   end
   
  % quiver(xs(1:elementSkip:idx-1), ys(1:elementSkip:idx-1), qux(1:elementSkip:idx-1), quy(1:elementSkip:idx-1), 4);  
   surf(hm, 'EdgeColor','none');

   az = 0;
el = 90;
view(az, el);
title('Function plot');
xlabel('x');
ylabel('y');

%figure
hold on 
axis([0 mat_cols 0 mat_cols])
drawnow
hold off
 end












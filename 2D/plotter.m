clear all;
pressurePlot = true;
mat_rows = 11*2*10 + 1;

if pressurePlot
    data = dlmread('pres');
else
    data = dlmread('out');
end

fan = dlmread('fan.txt');
fanN =  size(fan)
n = size(data);
 data_rows = n(1);
 data_cols = n(2);
 mat_cols = data_cols;
 skip = 1;
 elementSkip = 1;
 
 sys_rows = data_cols;
 sys_cols = data_cols;
 xc = sys_cols/2;
 yc = sys_rows/2; %revisar

  
fanIndex = 1
pause(2);
 for base=1:skip*mat_rows*2:data_rows
     U = data(base:base+mat_rows-1, 1:mat_cols);
     V = data(base+mat_rows: base+2*mat_rows-1, 1:mat_cols);
     base
     fanIndex = fanIndex+1;
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
               if pressurePlot
                    hm(j,i) = qux(idx);
               else
                    hm(j,i) = sqrt(qux(idx)*qux(idx) + quy(idx)*quy(idx));
               end
 
               idx = idx+1;
       end
   end
   
   %quiver(xs(1:elementSkip:idx-1), ys(1:elementSkip:idx-1), qux(1:elementSkip:idx-1), quy(1:elementSkip:idx-1), 1);  
   surf(hm, 'EdgeColor','none');

   az = 0;
el = 90;
view(az, el);
title('Function plot');
xlabel('x');
ylabel('y');
%pause(0.5)
%figure
hold on 
a0 = fan(fanIndex*skip, 3)
b0 = fan(fanIndex*skip, 4)
a1 = fan(fanIndex*skip, 1)
b1 = fan(fanIndex*skip, 2)

a0n = 2*yc-a0;
a1n = 2*yc-a1;
b0n = 2*yc-b0;
b1n = 2*yc-b1;

A = [a0n a1n] 
B = [b0n b1n] 
plot(A,B,'LineWidth',5)
axis([0 mat_cols 0 mat_cols])
hold on
line(A,B)
hold off
drawnow
hold off
 end












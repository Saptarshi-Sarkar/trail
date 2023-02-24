% Script to plot results
clc; clear;

D_c_array = 0.5:.5:2;
D_h_array = 0.2:.2:0.9;

[x, y] = ndgrid(D_c_array, D_h_array);

for i = 1:length(D_c_array)
    for j = 1:length(D_h_array)
        
        output = Nonlinear_analysis_func([D_c_array(i), D_h_array(j)]);
        
        wr(i,j) = output(1);
        zd(i,j) = output(2);
        fval(i,j) = output(3);
        beta(i,j) = output(4);
    end
end

figure(1)
surf(x,y,wr)
xlabel('D_c')
ylabel('D_h');
zlabel('\omega_r');

figure(2)
surf(x,y,zd)
xlabel('D_c')
ylabel('D_h');
zlabel('\zeta_d');

figure(3)
surf(x,y,fval)
xlabel('D_c')
ylabel('D_h');
zlabel('fval');

figure(4)
surf(x,y,beta)
xlabel('D_c')
ylabel('D_h');
zlabel('\beta');
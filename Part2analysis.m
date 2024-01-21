% Import data
time = Data2(:,1);
time = table2array(time);
%time = abs(time);
I = Data2(:,2);
I = table2array(I);
%I = abs(I);
Vx = Data2(:,3);
Vx = table2array(Vx);
%Vx = abs(Vx);
Vh = Data2(:,4);
Vh = table2array(Vh);
%Vh = abs(Vh);
Vh = 0.001*Vh;

% Standard Values
a = 0.0159;
b = 0.0095;
c = 0.00089;
err_abc = 0.0002;
q = 1.6E-19; 
I_high_err = 5.5E-03;
V_err = 3.1E-03;
Ix = 0.350;
err_Ix = 0.002;
u = 5000;
err_u = 0.02*u;
Rh = 7.12E-04;
err_Rh = 0.14E-04;


%Fit for ratio of N/L
H = [0.3,0.4,0.5];
H = transpose(H);
Ic = [0.899,1.369,2.028];
Ic = transpose(Ic);
err_H = 0.025*H;
err_Ic = 0.025*Ic;

[Ic_H,goF_Ic_H,Fit_output_Ic_H] = fit(Ic,H,'(m*x)+b','Weight',err_H.^(-2));
N_L = Ic_H.m;

    % Extract weighted jacobian
    J_Ic_H = Fit_output_Ic_H.Jacobian;
    
    %Get the covariance and curvature matrix and extract the errors on F
    %parameters from there.
    curvature_matrix_Ic_H = J_Ic_H'*J_Ic_H;
    covariance_matrix_Ic_H = inv(curvature_matrix_Ic_H);
    
    % Calculate CHI_squared
    min_chi2_Ic_H = goF_Ic_H.sse;
    dof_Ic_H = goF_Ic_H.dfe;
    reduced_chi2_Ic_H = min_chi2_Ic_H/dof_Ic_H;

    err_N_L = covariance_matrix_Ic_H(1,1);

%Calculate Bz = Vh*(c*Rh/Ix)

Bz = Vh.*(c*Rh/Ix);
err_Bz = (((Rh*c*V_err)/Ix).^2) + (((Vh*c*err_Rh)/Ix).^2) + (((Vh*Rh*c*err_Ix)/(Ix^2)).^2) + (((Vh*Rh*err_abc)/Ix).^2);
err_Bz = sqrt(err_Bz);

%Calculating H = (N_L)*I;
H = (u*I*(N_L)) + 60;
err_H = sqrt( ((N_L*u*I_high_err).^2) + ((I*u*err_N_L).^2) + ((I*N_L*err_u).^2) );

plot (H, Bz)
hold on
 xlabel('Magnetic Field Strength at the center of the Iron Core (A/m)','FontSize',18);
 ylabel('Magnetic Flux Density (T)','FontSize',18);
 grid on
 title ('Hysterisis Curve: Variation in Magnetic Flux Density and Magnetic Field Strength Based on direction of current ','FontSize',18)
 errorbar(H ,Bz, err_Bz,"k.",'HandleVisibility','off');
 errorbar(H ,Bz, err_H,'horizontal',"k.",'HandleVisibility','off');
 legend('Coresive Force is 60A/m, Remanent Flux Density is (1.42E-08)T','FontSize',18)
 hold off







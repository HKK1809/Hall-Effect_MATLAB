                   % Import Values
I3 = Data1(:,1);
I3 = table2array(I3);
I3 = abs(I3);
I4 = Data1(:,5);
I4 = table2array(I4);
I4 = abs(I4);
I5 = Data1(:,9);
I5 = table2array(I5);
I5 = abs(I5);

Vx3 = Data1(:,2);
Vx3 = table2array(Vx3);
Vx3 = abs(Vx3);
Vx4 = Data1(:,6);
Vx4 = table2array(Vx4);
Vx4 = abs(Vx4);
Vx5 = Data1(:,10);
Vx5 = table2array(Vx5);
Vx5 = abs(Vx5);

Vh3 = Data1(:,3);
Vh3 = table2array(Vh3);
Vh3 = 0.001*Vh3;
Vh3 = abs(Vh3);
Vh4 = Data1(:,7);
Vh4 = table2array(Vh4);
Vh4 = 0.001*Vh4;
Vh4 = abs(Vh4);
Vh5 = Data1(:,11);
Vh5 = table2array(Vh5);
Vh5 = 0.001*Vh5;
Vh5 = abs(Vh5);

a = 0.0159;
err_ab = 0.0002;
b = 0.0095;
c = 0.00089;
err_c = 0.00002;
q = 1.6E-19; 
err_T = 0.2;

err_B3 = 2.5*0.01*0.3;
err_B4 = 2.5*0.01*0.4;
err_B5 = 2.5*0.01*0.5;

 %Error
 I_high_err = (5.5E-03)*ones(size(I3));
 I_err = (0.31E-03)*ones(size(I3));
 V_err = (3.1E-03)*ones(size(I3));

 %Curve Fit

 %I Vs VH
 figure(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I3 Vs VH3
[I_Vh3,goF_I_Vh3,Fit_output_I_Vh3] = fit(I3,Vh3,'(m*x)+b','Weight',V_err.^(-2));

    % Extract weighted jacobian
    J_I_Vh3 = Fit_output_I_Vh3.Jacobian;
    
    %Get the covariance and curvature matrix and extract the errors on F
    %parameters from there.
    curvature_matrix_I_Vh3 = J_I_Vh3'*J_I_Vh3;
    covariance_matrix_I_Vh3 = inv(curvature_matrix_I_Vh3);
    
    % Calculate CHI_squared
    min_chi2_I_Vh3 = goF_I_Vh3.sse;
    dof_I_Vh3 = goF_I_Vh3.dfe;
    
    reduced_chi2_I_Vh3 = min_chi2_I_Vh3/dof_I_Vh3;
    
    err_I_Vh3_m = covariance_matrix_I_Vh3(1,1);

Rh_I_Vh3_m = (I_Vh3.m)*c/0.3;
err_Rh_I_Vh3_m = sqrt((((c*err_I_Vh3_m)/0.3)^2) +(((I_Vh3.m*err_c)/0.3)^2) + (((I_Vh3.m*c*err_B3)/0.09)^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I4 Vs VH4
[I_Vh4,goF_I_Vh4,Fit_output_I_Vh4] = fit(I4,Vh4,'(m*x)+b','Weight',V_err.^(-2));

    % Extract weighted jacobian
    J_I_Vh4 = Fit_output_I_Vh4.Jacobian;
    
    %Get the covariance and curvature matrix and extract the errors on F
    %parameters from there.
    curvature_matrix_I_Vh4 = J_I_Vh4'*J_I_Vh4;
    covariance_matrix_I_Vh4 = inv(curvature_matrix_I_Vh4);
    
    % Calculate CHI_squared
    min_chi2_I_Vh4 = goF_I_Vh4.sse;
    dof_I_Vh4 = goF_I_Vh4.dfe;
    
    reduced_chi2_I_Vh4 = min_chi2_I_Vh4/dof_I_Vh4;
    
    err_I_Vh4_m = covariance_matrix_I_Vh4(1,1);

Rh_I_Vh4_m = (I_Vh4.m)*c/0.4;
err_Rh_I_Vh4_m = sqrt((((c*err_I_Vh4_m)/0.4)^2) +(((I_Vh4.m*err_c)/0.4)^2) + (((I_Vh4.m*c*err_B4)/0.16)^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I5 Vs Vh5
[I_Vh5,goF_I_Vh5,Fit_output_I_Vh5] = fit(I5,Vh5,'(m*x)+b','Weight',V_err.^(-2));

    % Extract weighted jacobian
    J_I_Vh5 = Fit_output_I_Vh5.Jacobian;
    
    %Get the covariance and curvature matrix and extract the errors on F
    %parameters from there.
    curvature_matrix_I_Vh5 = J_I_Vh5'*J_I_Vh5;
    covariance_matrix_I_Vh5 = inv(curvature_matrix_I_Vh5);
    
    % Calculate CHI_squared
    min_chi2_I_Vh5 = goF_I_Vh5.sse;
    dof_I_Vh5 = goF_I_Vh5.dfe;
    
    reduced_chi2_I_Vh5 = min_chi2_I_Vh5/dof_I_Vh5;
    
    err_I_Vh5_m = covariance_matrix_I_Vh5(1,1);

Rh_I_Vh5_m = (I_Vh5.m)*c/0.5;
err_Rh_I_Vh5_m = sqrt((((c*err_I_Vh5_m)/0.5)^2) +(((I_Vh5.m*err_c)/0.5)^2) + (((I_Vh5.m*c*err_B5)/0.25)^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate value of Hall Constant and mean error on it
Rh = (Rh_I_Vh3_m + Rh_I_Vh4_m + Rh_I_Vh5_m)/3;
err_Rh = (sqrt((err_Rh_I_Vh3_m*err_Rh_I_Vh3_m)+(err_Rh_I_Vh4_m*err_Rh_I_Vh4_m)+(err_Rh_I_Vh5_m*err_Rh_I_Vh5_m)))/3;
 
n = 1/(q*Rh);
err_n = err_Rh/(q*Rh*Rh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLotting Figure 1
 F1_3 = plot(I_Vh3,'r',I3,Vh3);hold on;
 F1_4 = plot(I_Vh4,'g',I4,Vh4);
 F1_5 = plot(I_Vh5,'b',I5,Vh5);
 
 F1 =[F1_3;F1_4;F1_5];
 xlabel('Applied Current- Ix (A)','FontSize',18);
 ylabel('Hall Potential-Vh (V)','FontSize',18);
 title ('Variation in Hall Potential as the applied current is Varied at a constant magnetic field','FontSize',18)
 
 legend('Data for 0.3 Tesla', 'Rh = (6.79 +/- 0.23)E-04','Data for 0.4 Tesla','Rh = (7.14 +/- 0.24)E-04','Data for 0.5 Tesla','Rh = (7.43 +/- 0.25)E-04','FontSize',18)
 errorbar(I3 ,Vh3, I_err,'horizontal',"k.",'HandleVisibility','off');
 errorbar(I3 ,Vh3, V_err,"k.",'HandleVisibility','off');
 errorbar(I4 ,Vh4, I_err,'horizontal',"k.",'HandleVisibility','off');
 errorbar(I4 ,Vh4, V_err,"k.",'HandleVisibility','off');
 errorbar(I5 ,Vh5, I_err,'horizontal',"k.",'HandleVisibility','off');
 errorbar(I5 ,Vh5, V_err,"k.",'HandleVisibility','off');
 hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Vx Vs Vh
figure(2);
%Vx3 Vs Vh3
[Vx_Vh3,goF_Vx_Vh3,Fit_output_Vx_Vh3] = fit(Vx3,Vh3,'(m*x)+b','Weight',V_err.^(-2));

mobility3 = ((Vx_Vh3.m)*a)/(0.3*b);
    % Extract weighted jacobian
    J_Vx_Vh3 = Fit_output_Vx_Vh3.Jacobian;
    
    %Get the covariance and curvature matrix and extract the errors on F
    %parameters from there.
    curvature_matrix_Vx_Vh3 = J_Vx_Vh3'*J_Vx_Vh3;
    covariance_matrix_Vx_Vh3 = inv(curvature_matrix_Vx_Vh3);
    
    % Calculate CHI_squared
    min_chi2_Vx_Vh3 = goF_Vx_Vh3.sse;
    dof_Vx_Vh3 = goF_Vx_Vh3.dfe;
    
    reduced_chi2_Vx_Vh3 = min_chi2_Vx_Vh3/dof_Vx_Vh3;

err_Vx_Vh3_m = covariance_matrix_Vx_Vh3(1,1);
err_mobility_Vx_Vh3 = sqrt((err_Vx_Vh3_m*(a/(0.3*b)))^2 + (((Vx_Vh3.m)*err_ab)/(0.3*b))^2 + (((Vx_Vh3.m)*a*err_ab)/(0.3*b*b))^2 +  (((Vx_Vh3.m)*a*err_T)/(0.09*b))^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vx4 Vs Vh4
[Vx_Vh4,goF_Vx_Vh4,Fit_output_Vx_Vh4] = fit(Vx4,Vh4,'(m*x)+b','Weight',V_err.^(-2));

mobility4 = ((Vx_Vh4.m)*a)/(0.4*b);
    % Extract weighted jacobian
    J_Vx_Vh4 = Fit_output_Vx_Vh4.Jacobian;
    
    %Get the covariance and curvature matrix and extract the errors on F
    %parameters from there.
    curvature_matrix_Vx_Vh4 = J_Vx_Vh4'*J_Vx_Vh4;
    covariance_matrix_Vx_Vh4 = inv(curvature_matrix_Vx_Vh4);
    
    % Calculate CHI_squared
    min_chi2_Vx_Vh4 = goF_Vx_Vh4.sse;
    dof_Vx_Vh4 = goF_Vx_Vh4.dfe;
    
    reduced_chi2_Vx_Vh4 = min_chi2_Vx_Vh4/dof_Vx_Vh4;

err_Vx_Vh4_m = covariance_matrix_Vx_Vh4(1,1);
err_mobility_Vx_Vh4 = sqrt((err_Vx_Vh4_m*(a/(0.4*b)))^2 + (((Vx_Vh4.m)*err_ab)/(0.4*b))^2 + (((Vx_Vh4.m)*a*err_ab)/(0.4*b*b))^2 +  (((Vx_Vh3.m)*a*err_T)/(0.16*b))^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vx5 Vs Vh5
[Vx_Vh5,goF_Vx_Vh5,Fit_output_Vx_Vh5] = fit(Vx5,Vh5,'(m*x)+b','Weight',V_err.^(-2));

mobility5 = ((Vx_Vh5.m)*a)/(0.5*b);
    % Extract weighted jacobian
    J_Vx_Vh5 = Fit_output_Vx_Vh5.Jacobian;
    
    %Get the covariance and curvature matrix and extract the errors on F
    %parameters from there.
    curvature_matrix_Vx_Vh5 = J_Vx_Vh5'*J_Vx_Vh5;
    covariance_matrix_Vx_Vh5 = inv(curvature_matrix_Vx_Vh5);
    
    % Calculate CHI_squared
    min_chi2_Vx_Vh5 = goF_Vx_Vh5.sse;
    dof_Vx_Vh5 = goF_Vx_Vh5.dfe;
    
    reduced_chi2_Vx_Vh5 = min_chi2_Vx_Vh5/dof_Vx_Vh5;

err_Vx_Vh5_m = covariance_matrix_Vx_Vh5(1,1);
err_mobility_Vx_Vh5 = sqrt((err_Vx_Vh5_m*(a/(0.5*b)))^2 + (((Vx_Vh5.m)*err_ab)/(0.5*b))^2 + (((Vx_Vh5.m)*a*err_ab)/(0.5*b*b))^2 +  (((Vx_Vh5.m)*a*err_T)/(0.25*b))^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate mean mobility and error on it

mobility = (mobility3 + mobility4 + mobility5)/3;
err_mobility = (sqrt((err_mobility_Vx_Vh3*err_mobility_Vx_Vh3)+(err_mobility_Vx_Vh4*err_mobility_Vx_Vh4)+(err_mobility_Vx_Vh5*err_mobility_Vx_Vh5)))/3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potting figure 2
 VV3 = plot(Vx_Vh3,'r',Vx3,Vh3);hold on;
 VV4 =  plot(Vx_Vh4,'g',Vx4,Vh4);
 VV5 = plot(Vx_Vh5,'b',Vx5,Vh5);
 
 VV =[VV3;VV4;VV5];
 
 xlabel('Applied Potential- Vx (V)','FontSize',18);
 ylabel('Hall Potential-Vh (V)','FontSize',18);
 title ('Variation in Hall Potential as the applied Potential is Varied at a constant magnetic field','FontSize',18)
 legend('Data for 0.3 Tesla', 'mobility = (0.38 +/- 0.25)','Data for 0.4 Tesla','mobility = (0.39 +/- 0.14)','Data for 0.5 Tesla','mobility = (0.40 +/- 0.16)','FontSize',18)
 
 errorbar(Vx5, Vh5, V_err,"k.",'HandleVisibility','off');
 errorbar(Vx4, Vh4, V_err,"k.",'HandleVisibility','off');
 errorbar(Vx3, Vh3, V_err,"k.",'HandleVisibility','off');
  errorbar(Vx5, Vh5, V_err,'horizontal',"k.",'HandleVisibility','off');
 errorbar(Vx4, Vh4, V_err,'horizontal',"k.",'HandleVisibility','off');
 errorbar(Vx3, Vh3, V_err,'horizontal',"k.",'HandleVisibility','off');
 hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %I Vs Vx
 figure(3);
% Making Fit
 [I_Vx3,goF_I_Vx3,Fit_output_I_Vx3] = fit(I3,Vx3,'(m*x)+b','Weight',V_err.^(-2));
 [I_Vx4,goF_I_Vx4,Fit_output_I_Vx4] = fit(I4,Vx4,'(m*x)+b','Weight',V_err.^(-2));
 [I_Vx5,goF_I_Vx5,Fit_output_I_Vx5] = fit(I5,Vx5,'(m*x)+b','Weight',V_err.^(-2));
 % Calculating Resistance
 R3 = I_Vx3.m;
 R4 = I_Vx4.m; 
 R5 = I_Vx5.m;
% Calculating Resistivity
 resistivity3 = (R3*b*c)/a;
 resistivity4 = (R4*b*c)/a; 
 resistivity5 = (R5*b*c)/a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % error propogation
    %0.3Tesla Error Prpogation
     % Extract weighted jacobian
    J_I_Vx3 = Fit_output_I_Vx3.Jacobian;
    
    %Get the covariance and curvature matrix and extract the errors on F
    %parameters from there.
    curvature_matrix_I_Vx3 = J_I_Vx3'*J_I_Vx3;
    covariance_matrix_I_Vx3 = inv(curvature_matrix_I_Vx3);
    
    % Calculate CHI_squared
    min_chi2_I_Vx3 = goF_I_Vx3.sse;
    dof_I_Vx3 = goF_I_Vx3.dfe;
    
    reduced_chi2_I_Vx3 = min_chi2_I_Vx3/dof_I_Vx3;

    err_R3 = covariance_matrix_I_Vx3(1,1);
    err_resistivity3 = sqrt((((err_R3*b*c)/a)^2) + (((R3*c*err_ab)/a)^2) + (((R3*b*err_ab)/a)^2) + (((R3*b*c*err_ab)/(a*a))^2));
    
    
    %0.4Tesla Error Prpogation
     % Extract weighted jacobian
    J_I_Vx4 = Fit_output_I_Vx4.Jacobian;
    
    %Get the covariance and curvature matrix and extract the errors on F
    %parameters from there.
    curvature_matrix_I_Vx4 = J_I_Vx4'*J_I_Vx4;
    covariance_matrix_I_Vx4 = inv(curvature_matrix_I_Vx4);
    
    % Calculate CHI_squared
    min_chi2_I_Vx4 = goF_I_Vx4.sse;
    dof_I_Vx4 = goF_I_Vx4.dfe;
    
    reduced_chi2_I_Vx4 = min_chi2_I_Vx4/dof_I_Vx4;

    err_R4 = covariance_matrix_I_Vx4(1,1);
    err_resistivity4 = sqrt(((err_R4*b*c)/a)^2 + ((R4*c*err_ab)/a)^2 + ((R4*b*err_ab)/a)^2 + ((R4*b*c*err_ab)/(a*a))^2);

    %0.5 Tesla Error Prpogation
     % Extract weighted jacobian
    J_I_Vx5 = Fit_output_I_Vx5.Jacobian;
    
    %Get the covariance and curvature matrix and extract the errors on F
    %parameters from there.
    curvature_matrix_I_Vx5 = J_I_Vx5'*J_I_Vx5;
    covariance_matrix_I_Vx5 = inv(curvature_matrix_I_Vx5);
    
    % Calculate CHI_squared
    min_chi2_I_Vx5 = goF_I_Vx5.sse;
    dof_I_Vx5 = goF_I_Vx5.dfe;
    
    reduced_chi2_I_Vx5 = min_chi2_I_Vx5/dof_I_Vx5;

    err_R5 = covariance_matrix_I_Vx5(1,1);
    err_resistivity5 = sqrt(((err_R5*b*c)/a)^2 + ((R5*c*err_ab)/a)^2 + ((R5*b*err_ab)/a)^2 + ((R5*b*c*err_ab)/(a*a))^2);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % calculating Net resistivity
 resistivity = (resistivity3 + resistivity4 + resistivity5)/3;
 err_resistivity = ((err_resistivity3)^2) + ((err_resistivity4)^2) + ((err_resistivity5)^2);
 err_resistivity = sqrt(err_resistivity);
 err_resistivity = err_resistivity/3;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Plot Figure 3
  F3_3 = plot(I_Vx3,'r',I3,Vx3); hold on;
  F3_4 = plot(I_Vx4,'g',I4,Vx4);
  F3_5 = plot(I_Vx5,'b',I5,Vx5);
 
 F3 =[F3_3;F3_4;F3_5];
 xlabel('Applied Current- Ix (A)','FontSize',18);
 ylabel('Applied Potential-Vx (V)','FontSize',18);
 title ('Variation in Potential as the applied current is Varied at a constant magnetic field','FontSize',18)
 
  legend('Data for 0.3 Tesla', 'Resistivity = (18.1 +/- 4.1)E-04','Data for 0.4 Tesla','Resistivity = (18.3 +/- 4.1)E-04','Data for 0.5 Tesla','Resistivity = (18.7 +/- 4.2)E-04','FontSize',18)
  errorbar(I3 ,Vx3, I_err,'horizontal',"k.",'HandleVisibility','off');
  errorbar(I3 ,Vx3, V_err,"k.",'HandleVisibility','off');
  errorbar(I4 ,Vx4, I_err,'horizontal',"k.",'HandleVisibility','off');
  errorbar(I4 ,Vx4, V_err,"k.",'HandleVisibility','off');
  errorbar(I5 ,Vx5, I_err,'horizontal',"k.",'HandleVisibility','off');
  errorbar(I5 ,Vx5, V_err,"k.",'HandleVisibility','off');
  hold off









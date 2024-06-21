
function [col_1, col_2, col_e, f, p] = rec_cylindrical(D, r_max, r_min, V, parameters, ind_screening, polarity, tipo_gamma)


%%% function to study recombination in cylindrical ionization chambers

%INPUTS
%D: dose [Gy]
%r_min: inner radius of the ionization chamber [m]
%r_max: outer radius of the ionization chamber [m]
%V: operation voltage [V]
%parameters: vector of physical parameters
%ind_screening: flag controlling the calculation or not of field screening effects; if ind_screening=1 we include field screening in the computation
%polarity: flag controlling the detector polarity; if polarity=1 negative charge moves towards r_max, if polarity=-1 negative charge moves towards r_min
%tipo gamma: flag controlling the calculation of gamma; if tipo_gamma=1, we include the electric field dependence in gamma, otherwise we consider gamma as a constant 


%OUTPUTS
%col_1: collection of positive ions
%col_2: collection of negative ions
%col_e: collection of electrons
%f: collection efficiency
%p: free electron fraction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INITIALIZATION

k1=parameters(1); %ion mobility [m^2 V^(-1) s^(-1)]
k2=parameters(2); %ion mobility [m^2 V^(-1) s^(-1)] 
alpha=parameters(3); %recombination constant [m^3 s^(-1)] 

ke=parameters(4); %electron mobility [m^2 V^(-1) s^(-1)]


%space discretization
dx=0.5e-6;

%spatial coordinates and electric field
x=[r_min+dx:dx:r_max] - dx/2;
E=V./(log(r_max/r_min)*x);

E0=E;

%time discretization
dt_i=dx/(25*k2*max(E));
dt_e=dx/(25*ke*max(E));

%initial clculation of gamma
switch tipo_gamma
    case 1
        gamma = gamma_dep_E(E);
    
    otherwise 
        valor=2e5;
        gamma = gamma_dep_E(valor);
end


ind=length(x);

%initial ionization density in air [pairs/m^3]
n0=2.216e17*D;
n1(1:ind)=n0;
n2(1:ind)=0;
ne(1:ind)=n0;



total_charge=2*sum(x.*n1);
total_rec=0;



t=0;

col_1=0;
col_2=0;
col_e=0;

%flag that controls the change of regime from "free electrons" to "no free electrons"
cambio=0;

dt=dt_e;



%DYNAMICS
while mean(n1)/n0>1e-6
    

    
    %recombination
    rec=alpha*n1.*n2*dt;
    %we ignore recombination of electrons
    %rec_e=alpha_e*n1.*ne*dt;
    %ne=ne-rec_e;
    

    %attachment
    att=gamma.*ne*dt;
    
    
    n1=n1-rec;
    n2=n2-rec+att;
    ne=ne-att;
    
    total_rec = total_rec + 2*sum(x.*rec);
   
    %avoid negative numbers that may arise for inadequate discretizations
    find(n1<0);, n1(ans)=0;
    find(n2<0);, n2(ans)=0;
    find(ne<0);, ne(ans)=0;
        
    n1_new=n1;
    n2_new=n2;
    ne_new=ne;
    

    %DRIFT
    delta_x1=k1*E*dt;
    delta_x2=k2*E*dt;
    delta_xe=ke*E*dt;
   
    for i=1:ind
       
        
        switch polarity
            
            case 1 
                %positive, positive charge moves towards r_min
        
                if i==1
                    
                    n1_new(i)=n1(i)*(1-delta_x1(i)/dx) + (x(i+1)/x(i))*n1(i+1)*(delta_x1(i+1)/dx);
                    n2_new(i)=n2(i)*(1-delta_x2(i)/dx);

                    ne_new(i)=ne(i)*(1-delta_xe(i)/dx);

                    col_1=col_1+n1(i)*(delta_x1(i)/dx)*x(i);

                elseif i==ind

                    n1_new(i)=n1(i)*(1-delta_x1(i)/dx);
                    n2_new(i)=n2(i)*(1-delta_x2(i)/dx) + (x(i-1)/x(i))*n2(i-1)*(delta_x2(i-1)/dx);

                    ne_new(i)=ne(i)*(1-delta_xe(i)/dx) + (x(i-1)/x(i))*ne(i-1)*(delta_xe(i-1)/dx);

                    col_2=col_2+n2(i)*(delta_x2(i)/dx)*x(i);
                    col_e=col_e+ne(i)*(delta_xe(i)/dx)*x(i);

                else

                    n1_new(i)=n1(i)*(1-delta_x1(i)/dx) + (x(i+1)/x(i))*n1(i+1)*(delta_x1(i+1)/dx);
                    n2_new(i)=n2(i)*(1-delta_x2(i)/dx) + (x(i-1)/x(i))*n2(i-1)*(delta_x2(i-1)/dx);

                    ne_new(i)=ne(i)*(1-delta_xe(i)/dx) + (x(i-1)/x(i))*ne(i-1)*(delta_xe(i-1)/dx);

                end
                
                
            case -1
                
                %negative, negative charge moves towards r_min
                
                if i==1
                    n2_new(i)=n2(i)*(1-delta_x2(i)/dx) + (x(i+1)/x(i))*n2(i+1)*(delta_x2(i+1)/dx);
                    n1_new(i)=n1(i)*(1-delta_x1(i)/dx);
                    
                    ne_new(i)=ne(i)*(1-delta_xe(i)/dx) + (x(i+1)/x(i))*ne(i+1)*(delta_xe(i+1)/dx);

                    %charge collection
                    col_2=col_2+n2(i)*(delta_x2(i)/dx)*x(i);
                    col_e=col_e+ne(i)*(delta_xe(i)/dx)*x(i);

                elseif i==ind

                    n2_new(i)=n2(i)*(1-delta_x2(i)/dx);
                    n1_new(i)=n1(i)*(1-delta_x1(i)/dx) + (x(i-1)/x(i))*n1(i-1)*(delta_x1(i-1)/dx);
                    
                    ne_new(i)=ne(i)*(1-delta_xe(i)/dx);

                    %charge collection
                    col_1=col_1+n1(i)*(delta_x1(i)/dx)*x(i);


                else
                    n2_new(i)=n2(i)*(1-delta_x2(i)/dx) + (x(i+1)/x(i))*n2(i+1)*(delta_x2(i+1)/dx);
                    n1_new(i)=n1(i)*(1-delta_x1(i)/dx) + (x(i-1)/x(i))*n1(i-1)*(delta_x1(i-1)/dx);
                    
                    ne_new(i)=ne(i)*(1-delta_xe(i)/dx) + (x(i+1)/x(i))*ne(i+1)*(delta_xe(i+1)/dx);

                    
                end
                
        end
                
        
    end


    n1=n1_new;
    n2=n2_new;
    ne=ne_new;

    %avoid negative numbers that may arise for inadequate discretizations
    find(n1<0);, n1(ans)=0;        
    find(n2<0);, n2(ans)=0;
    find(ne<0);, ne(ans)=0;
        
       

    %field screening
    if ind_screening==1

        E_ap(1:ind)=0;
        for j=1:ind-1
            E_ap(j+1)=E_ap(j)-dx*((1.602e-19/(1.0006*8.854e-12))*polarity*(n1(j)-n2(j)-ne(j)));
        end
        E=E_ap+E0;

        %effect of power supply
        E = E + (V-sum(E)*dx)./(sum(1./x)*x*dx);


        find(E<0);, E(ans)=0;

        %recomputation of gamma due to field screening
        if tipo_gamma==1
            gamma = gamma_dep_E(E);
        end

        if cambio==1
            dt=dx/(25*k2*max(E));
        elseif cambio==0
            dt=dx/(25*ke*max(E));
        end

    end
 
    t=t+dt;
    
    %if "almost all" electrons have been collected we update the value of dt
    if mean(ne)/n0<1e-9 && cambio==0
        dt=dt_i;
        cambio=1;
        col_e=col_e+sum(ne.*x);
        ne=ne*0;
    end
    
end


col_1=col_1+sum(n1.*x);
col_2=col_2+sum(n2.*x);
col_e=col_e+sum(ne.*x);

%collection efficiency, defined as the ratio of collected change to produced charge
f=(col_1+col_2+col_e)/(total_charge);

%free electron fraction, defined as the ratio of collected electrons to ionized electrons
p=col_e/(total_charge/2);

%NOTE: formally, the readout signal is given not by the collected charge but by the induced charge in the
%electrodes (the Shockley-Ramo theorem, see Paz-Martin et al Phys Med 2022). For the purpose of this
%work, the approximation of using the collected charge is good enough.


end



function gamma = gamma_dep_E(E)

ind=find(E<0.35e5);
gamma=(1.1+11.3*exp(-1.04e-5*E))*1e7;
gamma(ind)=7e7+657*E(ind);

end



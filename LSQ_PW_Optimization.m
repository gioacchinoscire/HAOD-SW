function [sol_struct]=LSQ_PW_Optimization(X0,kmax,obs_pos,Cel_Coord,time,order,EGM,EOP,DAT,C_r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves the Nonlinear LSQ problem with the Powell's dog-leg 
% algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% X0: Initial Guess-Dimensional
% kmax: max iteration tolerance
% Cel_Coord: measurements matrix mx2
% time: siGMlation time JD
% C_lm,S_lm: Harmonics coefficients for the gravitational model EGM2008
% OUTPUT:
% Sol_struct. structure that contrains the solution data (optiGMm state,
% cost function, covariance matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global R_E GM

adim=[(R_E/1E3).*ones(3,1);sqrt((GM/1E9)/(R_E/1E3))*ones(3,1)];

Cost_fun=zeros(kmax,1);
Covariance=zeros(6,6,kmax);

eps1=1.e-16;
eps2=1.e-16;
eps3=1.e-16;


k=0; x=X0./adim;

fx=Residual_vec(x.*adim,time,obs_pos,order,EGM,EOP,DAT,Cel_Coord,C_r);
J=Ad_Jacobian(x,obs_pos,Cel_Coord,time,order,EGM,EOP,DAT,C_r);
g=J'*fx;
DELTA=1E0;

found=(norm(g,inf)<=eps1);

fprintf('Iteration \tobjective_function\n\n');

while found==0 && k<kmax
    k=k+1;
    alpha=(norm(g)/(norm(J*g)))^2;
    h_sd=-g;
    h_gn=(J'*J)\h_sd;
    Cost_fun(k)=ObjFunc(fx);
    fprintf('\t%i\t\t\t%d\n',k,Cost_fun(k));
    %% h_dl step determination.............................................
    
    if norm(h_gn)<=DELTA
        h_dl=h_gn;
        D_L_model=Cost_fun(k);
    elseif norm(alpha.*h_sd)>=DELTA
        h_dl=(DELTA/norm(h_sd)).*h_sd;
        D_L_model=DELTA*(2*norm(alpha.*g)-DELTA)/(2*alpha);
    else
        a=alpha.*h_sd;
        b=h_gn;
        c=a'*(b-a);
        if c<=0
            beta=(-c+sqrt(c^2+(norm(b-a)^2)*(DELTA^2-norm(a)^2)))/...
                (norm(b-a)^2);
        else
            beta=(DELTA^2-norm(a)^2)/(c+sqrt(c^2+(norm(b-a)^2)*...
                (DELTA^2-norm(a)^2)));
        end
        h_dl=alpha.*h_sd+beta*(h_gn-alpha.*h_sd);
        D_L_model=(1/2)*alpha*((1-beta)^2)*norm(g)^2+beta*(2-beta)*...
            Cost_fun(k);
    end
    %% minimization process................................................
    
    if norm(h_dl)<=eps2*(norm(x)+eps2)
        found=1;
    else
        x_new=(x+h_dl(1:6));
        C_r_new=C_r+h_dl(7);
        fx_new=Residual_vec(x_new.*adim,time,obs_pos,order,EGM,EOP,DAT,Cel_Coord,C_r_new);
        rho=(Cost_fun(k)-ObjFunc(fx_new))/D_L_model;
        if rho>0
            x=x_new;
            C_r=C_r_new;
            fx=fx_new;
            J=Ad_Jacobian(x,obs_pos,Cel_Coord,time,order,EGM,EOP,DAT,C_r);
            g=J'*fx;
            found=((norm(g,inf)<=eps1)||(norm(fx,inf)<=eps3));
            if rho>0.75
                DELTA=max([DELTA,3*norm(h_dl)]);
            end
        else
            DELTA=DELTA/2;
        end
    end   
    

    
end

X_opt=x.*adim;

Cost_fun=nonzeros(Cost_fun);

if k<kmax
    Covariance=Covariance(:,:,1:k+1);
else
    Covariance=Covariance(:,:,1:k);
end
    
H=observ_grad(X_opt,time,obs_pos,order,EGM,EOP,DAT,Cel_Coord,C_r);
R=(3/3600)^2*eye(numel(Cel_Coord));
Covariance(:,:,end)=inv(H'*(inv(R))*H);

sol_struct=struct('sol',X_opt,'Cr',C_r,'costfun',Cost_fun,'cov',Covariance);

end


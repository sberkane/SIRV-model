function [S,I,R,V] = SIRV_model(alpha,beta,gamma,I0,R0,V0,T,dt,u)
    S = zeros(1,ceil(T/dt));
    R = zeros(1,ceil(T/dt));
    I = zeros(1,ceil(T/dt));
    V = zeros(1,ceil(T/dt));
    S(1) =1-I0-R0-V0; 
    I(1) = I0;
    R(1)=R0;
    V(1)=V0;
    for t = 1:ceil(T/dt)-1
        %Boundary control
        if S(t)<=0
            u(t)=0;
            S(t)=0;
        end
        % Equations of the model
        dS = (-beta*I(t)*S(t)-alpha*u(t)) * dt;
        dI = (beta*I(t)*S(t) - gamma*I(t)) * dt;
        dR = (gamma*I(t)) * dt;
        dV=alpha*u(t) * dt;
        S(t+1) = S(t) + dS;
        I(t+1) = I(t) + dI;
        R(t+1) = R(t) + dR;
        V(t+1) = V(t) + dV;

    end
end
function [Svec1,Ivec1,Rvec1,x,timeStep,timeElapsed]=Euler_SIR(S0,I0,R0,a,r,days,h,tol)
    % Checks time elapsed
    tic;
    
    % Time step counter
    i = 1;
    
    % Set errors to infinity so that no matter the tolerance level, the
    % while loop will always execute
    sError = inf;
    iError = inf;
    rError = inf;
    
    % While loop executes until the error in susceptibles, infectives, and
    % removed populations falls below tolerance level set by user
    while (sError>tol) || (iError>tol) || (rError>tol)
        % Reset result vectors at the start of every iteration
        Svec1=[S0];
        Ivec1=[I0];
        Rvec1=[R0];
        Svec2=[S0];
    	Ivec2=[I0];
        Rvec2=[R0];
        % Declare initial conditions
        S = S0;
        I = I0;
        R = R0;
        % First loop uses step size h/i where i is the time step counter
        for j=1:(days-1)*i
            Snew = S-(h/i)*r*I*S;
            Inew = I+(h/i)*(r*I*S-a*I);
            Rnew = R+(h/i)*(a*I);
            Svec1 = [Svec1; Snew];
            Ivec1 = [Ivec1; Inew];
            Rvec1 = [Rvec1; Rnew];
            S = Snew;
            I = Inew;
            R = Rnew;
        end
        % Increment time step counter
        i = i+1;
        % Reset initial conditions
        S = S0;
        I = I0;
        R = R0;
        % Second loop uses step size h/(i+1) to compare results from
        % previous for loop to results with results of higher accuracy
        for k=1:(days-1)*i
           Snew = S-(h/i)*r*I*S;
           Inew = I+(h/i)*(r*I*S-a*I);
           Rnew = R+(h/i)*(a*I);
           Svec2 = [Svec2; Snew];
           Ivec2 = [Ivec2; Inew];
           Rvec2 = [Rvec2; Rnew];
           S = Snew;
           I = Inew;
           R = Rnew; 
        end
        
        % Forward errors defined as the difference between the final values
        % of each population using time step h/i and using time step
        % h/(i+1) to provide operational definition for termination
        % condition
        sError = abs(Svec2(length(Svec2),1)-Svec1(length(Svec1),1));
        iError = abs(Ivec2(length(Ivec2),1)-Ivec1(length(Ivec1),1));
        rError = abs(Rvec2(length(Rvec2),1)-Rvec1(length(Rvec1),1));
        
    end
    
    timeElapsed = toc;
    % Return timeStep (must divide x-axis values by this number to obtain
    % time in days)
    timeStep = i-1;
    
    % Plotting results
    x=linspace(1,length(Svec1),length(Svec1))/timeStep;
    %plot(x,Svec1);
    %hold on
    %plot(x,Ivec1);
    %plot(x,Rvec1);
    %hold off
end
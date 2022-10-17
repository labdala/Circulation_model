function [t,y] = backward_euler_vec(f, tspan, ic , nsteps)
n = length(tspan);
h=(tspan(n)-tspan(1))/nsteps;
y=zeros(nsteps+1,length(ic));
t=zeros(nsteps+1,1);
y(1,:)=ic;
t(1)=tspan(1);
options = optimset('disp','off');
for i=2:(nsteps+1)
    t(i)=t(i-1)+h;
    y(i,:)=fsolve(@(Y) y(i-1,:)+h*f(t(i), Y)' - Y, y(i-1,:), options);
end

end
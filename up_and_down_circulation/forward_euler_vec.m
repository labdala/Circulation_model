function [t,y] = forward_euler_vec(f, tspan, ic , nsteps)
n = length(tspan);
h=(tspan(n)-tspan(1))/nsteps;
y=zeros(nsteps+1,length(ic));
t=zeros(nsteps+1,1);
y(1,:)=ic;
t(1)=tspan(1);
for i=2:(nsteps+1)
    t(i)=t(i-1)+h;
    %yi = y(i-1,:);
    %fi = f(t(i), y(i-1,:));
    y(i,:)=y(i-1,:)+h*f(t(i-1), y(i-1,:))';
end

end
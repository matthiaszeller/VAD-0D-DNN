function m = integralmean(y,t)
m=sum((y(1:end-1)+y(2:end)).*(t(2:end)-t(1:end-1))/2)/(t(end)-t(1));
end
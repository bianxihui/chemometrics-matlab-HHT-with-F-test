function [fnormhat,t]=instfreq(x,t,L,trace)
if (nargin == 0), 
error('At least one parameter required');
end; 
[xrow,xcol] = size(x); 
if (xcol~=1), 
error('X must have only one column'); 
end
if (nargin == 1),
    t=2:xrow-1;L=1;trace=0.0; 
elseif (nargin == 2), 
    L = 1;trace=0.0; 
elseif (nargin == 3), 
    trace=0.0; 
end;
if  L<1, 
    error('L must be >=1');
end 
[trow,tcol] = size(t);
if (trow~=1), 
error('T must have only one row'); 
end; 
if (L==1), 
  if any(t==1)||any(t==xrow), 
  error('T can not be equal to 1 neither to the last element of X');
  else 
  fnormhat=(pi-angle(-x(t).*conj(x(t-1))))/(2*pi); 
  end; 
else 
  H=kaytth(L); 
  if any(t<=L)||any(t+L>xrow), 
  error('The relation L<T<=length(X)-L must be satisfied'); 
  else
      for icol=1:tcol,
      if trace, disprog(icol,tcol,10); end; 
      ti = t(icol); tau = 0:L; 
      R = x(ti+tau).*conj(x(ti-tau));
      M4 = R(2:L+1).*conj(R(1:L));
      diff=2e-6;
      tetapred = H * (unwrap(angle(-M4))+pi); 
      while tetapred<0.0 , tetapred=tetapred+(2*pi);end;
      while tetapred>2*pi, tetapred=tetapred-(2*pi);end; 
      iter = 1; 
          while(diff > 1e-6)&&(iter<50),
             M4bis=M4.* exp(-j*2.0*tetapred); 
             teta = H * (unwrap(angle(M4bis))+2.0*tetapred); 
          while teta<0.0 , teta=(2*pi)+teta;end; 
          while teta>2*pi, teta=teta-(2*pi);end;   
          diff=abs(teta-tetapred); 
          tetapred=teta;
          iter=iter+1; 
          end; 
       fnormhat(icol,1)=teta/(2*pi); 
       end; 
   end; 
end;
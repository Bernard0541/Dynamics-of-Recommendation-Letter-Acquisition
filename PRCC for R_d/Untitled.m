%%%%%%%% Forward bifurcation parameters with Ro %%%%%%%%%%

   Pi=1.5;
   mu=0.01866;
   theta=0.854 ;
   rho1=0.805 ;
   eps= 0.3 ;
   gamma= 0.351 ;
   rho2= 0.219 ;
   m= 0.105; %0 .169795;
   delta=0.156 ;
   rho3=0.288;
   
  %mstar= (Lambda*alpha1)/(omega*mu);
  Q1= theta + rho1 + mu;
  Q2= mu + gamma + rho2;
  Q3= rho3+ mu;
  Q3= delta+ mu;
 % Ro= beta*c*gamma*(p*theta + (1-p)*Q2)*(1-mstar)./ Q1*Q2*Q3
  r=linspace(0,2,100);
  psy1=(eps*theta)/Q1;
  psy2=((1-eps)*theta + (gamma* psy1))/Q3;
  psy3=(rho1 + rho2*psy1 + rho3*psy2)/Q3;
  phi1= 1 + psy1 + psy2 + psy3 ;
  %phi2= beta*(1 + eta1*psy1 + eta2*psy2 + eta3*psy3) ;
  r_lst=linspace(0,2,100);
  mylist1=zeros(1,100);
  mylist2=zeros(1,100);
  ksi= Pi*Q1;
  %Ro= beta*c*gamma*(p*theta + (1-p)*Q2)*(1-mstar)./ Q1*Q2*Q3
  % B=omega*mu*Q1*Q2*Q3*(1-Ro);
  % root1=(-A1 + (A1^2 - 4*A2.*B).^(1/2)) ./ A2
  % root2=(-A1 - (A1^2 - 4*A2.*B).^(1/2)) ./ A2
for i=1:1:100
    r=r_lst(i);
       A2= Pi*Q1*(1-r);
       A1= Q1^2*(1-r) - Q1*phi1*(m*Pi + mu);
      Ao= m*Q1*phi1*(Q1 - mu*phi1);
    Rc=1 - (A1^2 ./(4*ksi.*Ao));
    mypoly=[Ao,A1,A2];
    sol=roots(mypoly);
    mylist1(i)=sol(1);
    mylist2(i)=sol(2);
    if mylist1(i)<0 || imag (mylist1(i))~=0
     mylist1(i)=0;
    end
    if mylist2(i)<0 || imag (mylist2(i))~=0
     mylist2(i)=0;
     
    end
end

Rc
mylist1;
mylist2;
k=1;
while mylist2(k)==0
    k=k+1;
end
ValD1=k;

m=ValD1;
while mylist2 (m)~=0
    m=m+1;
end
ValD2=m-1;
zer=zeros(1,ValD2);
s= linspace(0,1,100);
z= linspace(1,2,100);
y1=+ 0.2 + 0*s;
y2=+ 0.2 + 0*z;
% axis([0,2,0.2,3])
hold on
plot(r_lst(ValD1:end),mylist1(ValD1:end),'r',r_lst(ValD1:end),mylist2(ValD1:end),'b--','linewidth',1)
% plot(r_lst(ValD1:end),mylist1(ValD1:end),'r','linewidth',3)
plot(r_lst(1:ValD2),zer,'r','linewidth',1)
% plot(s, y1,'r','linewidth',2)
% plot(z, y2,'b--','linewidth',2)
xlabel('R_0')
ylabel('Infected Population')
 hold off
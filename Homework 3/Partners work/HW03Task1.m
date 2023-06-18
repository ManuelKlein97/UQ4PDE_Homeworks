%{
 I=4;
 h=1/I;
 nu=0.5;
 M0=10000;%100000;
 L=3;
 C=zeros(L+1,1);
 V=zeros(L+1,1);
 
 [C(1) V(1)]=GetConstants(@(x) CalculateQoI(0,L,x,I,h,nu),78,M0)
 for i=1:L
     i
     f=@(x)(CalculateQoI(i,L,x,2*I,1/(2*I),nu)-CalculateQoI(i,L,x,I,1/I,nu));
     [C(i+1) V(i+1)]=GetConstants(@(x) f(x),78,M0);
     I=2*I;
 end
 I=2;
 h=1/I;
 C
 V
 sum=0;
 for l=0:L
     sum=sum+sqrt(C(l+1)*V(l+1));
 end
 M=ones(3,L+1)*sum;
 M(1,:)=M(1,:)*0.1^(-2);
 M(2,:)=M(2,:)*0.05^(-2);
 M(3,:)=M(3,:)*0.01^(-2);
 for l=1:L
     M(:,l)=M(:,l)*sqrt(V(l)/C(l));
 end
 M=ceil(M)

M(1,3)

mlmc(@mlmc_l,M(1,L+1),0.1,L,10,0,0,0)
%}
 res=zeros(3,1);
 for i=1:3
     res(i,1)=MLMC(@(h,y) CalculateQoI(y,1/h,h,nu),M(i,:),h,2,L,78);
 end
 res
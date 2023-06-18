%{ 
I=4;
 h=1/I;
 nu=0.5;
 M0=10000;%100000;
 L=3;
 C=zeros(L,1);
 V=zeros(L,1);
 
 [C(1) V(1)]=GetConstants(@(x) CalculateQoI(x,I,h,nu),78,M0);
 for i=2:L
     i
     f=@(x)(CalculateQoI(x,2*I,1/(2*I),nu)-CalculateQoI(x,I,1/I,nu));
     [C(i) V(i)]=GetConstants(@(x) f(x),78,M0);
     I=2*I;
 end
 I=10;
 h=1/I;
 C
 V
 sum=0;
 for l=1:L
     sum=sum+sqrt(C(l)*V(l));
 end
 M=ones(3,L)*sum;
 M(1,:)=M(1,:)*0.1^(-2);
 M(2,:)=M(2,:)*0.05^(-2);
 M(3,:)=M(3,:)*0.01^(-2);
 for l=1:L
     M(:,l)=M(:,l)*sqrt(V(l)/C(l));
 end
 M=ceil(M)
%}
mlmc(@mlmc_l,M(1,1),0.1,2,4,0,0,0)

 res=zeros(3,1);
 for i=1:3
     res(i,1)=MLMC(@(h,y) CalculateQoI(y,1/h,h,nu),M(i,:),h,2,L,78);
 end
 res
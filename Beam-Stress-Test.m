clc
clear all

E0=200*10^9; I0=29*10^-6; h=0.15; l=1; F0=60*10^3;
q0=2.4*10^3;        % on introduit q0 dans le cas de q(x)=q0*x/l (cas du tp) 
nh=2;               % sinon on introduit le vecteur q (de meme dimension que le maillage)
                    % de meme pour E et I
function y=fpsi(i,x,he)
  if i==1
    y=1-3*x.^2/he^2 +2*x.^3/he^3;
  elseif i==2
    y=x-2*x.^2/he+x.^3/he^2;
  elseif i==3
    y=3*x.^2/he^2-2*x.^3/he^3;
  elseif i==4  
    y=-x.^2/he+x.^3/he^2;
  endif
endfunction


function y=psi2d(i,x,he)  % 2nd derivative
  if i==1
    y=-6/he^2+12/he^3*x;
  elseif i==2
    y=-4/he+6/he^2*x;
  elseif i==3
    y=6/he^2-12/he^3*x;
  elseif i==4
    y=-2/he+6*x/he^2;
  endif
endfunction


function [k]=kele(x1,x2,E,I)
  he=x2-x1;
  k=zeros(4,4);
  x=linspace(0,he,101);    
  for i=1:4
    for j=1:4
      k(i,j)=trapz( x, E.*I.*psi2d(i,x,he).*psi2d(j,x,he) );
    endfor
  endfor
endfunction

function [f]=fele(x1,x2,q) 
  he=x2-x1;
  f=zeros(1,4)';
  x=linspace(0,he,101);    
  for i=1:4
    f(i)=trapz( x, q.*fpsi(i,x,he) );
  endfor
endfunction


function x=maillage_regulier(nh,l)
  for i=1:nh+1
    x(i)=(i-1)*l/nh;
  endfor
endfunction


function x=maillage_variable(nh,l)
  for i=1:nh+1
    x(i)=l/2*(1-cos((i-1)/nh*pi));
  endfor
endfunction


function [k,f]=assemblage(nh,l,E,I,q,m)
  if m==0
    x=maillage_regulier(nh,l); 
  elseif m==1
    x=maillage_variable(nh,l); 
  endif
  k=zeros(2*nh+2);
  m=0;n=1;
  p=0;
  while p<2*nh
    m=m+1; n=n+1;
    x1=x(m);x2=x(n);
    z=101;             % on prend la partie de E et I qui appartient a l'element
    Eele=linspace(E(m),E(n),z); 
    Iele=linspace(I(m),I(n),z); 
    for i=1:4
      for j=1:4
        k(p+i,p+j)=k(p+i,p+j)+kele(x1,x2,Eele,Iele)(i,j);
      endfor
    endfor
    p=p+2;
  endwhile
  % de meme pour f
  f=zeros(2*nh+2,1);
  m=0;n=1;
  p=0;
  while p<2*nh
    m=m+1; n=n+1;
    x1=x(m);x2=x(n);
    z=101;                         % on prend la partie de q qui appartient a l'element
    qele=linspace(q(m),q(n),z);    % z represente le nombre d'element du vecteur q qui doit etre le meme que celui de x dans fele()
    for i=1:4
      f(p+i)=f(p+i)+fele(x1,x2,qele)(i);
    endfor
    p=p+2;
  endwhile
endfunction


function U=dep(nh,l)
  U=zeros(2*nh+2,1);
  U(1)=0;      % encastrement en x=0
  U(2)=0;
  for i=3:2*nh+2
    U(i)=NaN; %  la valeur NaN servira a reconnaitre les valeurs inconnues lors de la resolution
  endfor
endfunction


function Q=efforts(nh,l,F0)        % cond aux lims de la var secondaire
  Q=zeros(2*nh+2,1);
  Q(1)=NaN;
  Q(2)=NaN;
  for i=3:2*nh
    Q(i)=0;
  endfor
  Q(2*nh+1,1)=F0;
  Q(2*nh+2,1)=0;
endfunction

#Q=efforts(nh,l,F0)

function [vect,mat,eff]=inconnu(U,K,F)         % retourne le vecteur compose uniquement
  n=length(U);                                 % des valeurs inconnues et la matrice k
  vect=U; mat=K; eff=F;                        % correspondante
  while n>0
    if isnan(vect(n))==0   % on enleve les valeurs != nan qui sont connues
      vect(n)=[];
      eff(n)=[];
      mat(:,n)=[];
      mat(n,:)=[];
    endif
    n=n-1;
  endwhile
endfunction


function U=resolutionU(U,K,F)
  [U0,K0,F0]=inconnu(U,K,F);
  U0=K0\F0;                 % on calcule les valeurs inconnues de U
  j=1;
  for i=1:length(U)         % on met les nouvelles valeurs calculees de U
    if isnan(U(i))==1       % dans les composantes NaN
      U(i)=U0(j);
      j=j+1;
    endif
  endfor
endfunction


function Q=Qele(e,u,x,E,I,q)            % retourne le vecteur Q de l'element e.
  x1=x(e);                              % on doit recalculer K et F de l'element e
  x2=x(e+1);                            % et utiliser la partie du deplacement u
  he=x2-x1;                             % relative a l'element en question.
  z=101;         
  Eele=linspace(E(e),E(e+1),z); % on prend la partie de E, I et q qui appartient a l'element
  Iele=linspace(I(e),I(e+1),z); 
  qele=linspace(q(e),q(e+1),z);
  k=kele(x1,x2,Eele,Iele);          
  f=fele(x1,x2,qele);
  uele=u(2*(e-1)+1:2*(e-1)+1+3);
  Q=k*uele-f;
endfunction


function v=vAnalytique(x,L,E,I,F0,q0)
  v=(F0/(6*E*I)*x.^2).*(3*L-x)+q0*L^4/(120*E*I)*(20*x.^2/L^2-10*x.^3/L^3+x.^5/L^5);
endfunction


function main()
  E0=200*10^9; I0=29*10^-6; h=0.15; l=1; F0=60*10^3; q0=2.4*10^3;
  nh=8
  
  % maillage variable
  x=maillage_variable(nh,l);
  q=x*q0/l;
  for i=1:length(x)
    E(i)=E0;
    I(i)=I0;
  endfor
  [k,f]=assemblage(nh,l,E,I,q,1);
  u=dep(nh,l);
  Q=efforts(nh,l,F0);
  u=resolutionU(u,k,Q+f);
  % pour dessiner le graphique on prend uniquement les valeurs de u relatives a des translations
  j=1;
  for i=1:2*nh+2
    if mod(i,2)==1
      Ut(j)=u(i);
      j=j+1;
    endif
  endfor
  v=vAnalytique(x,l,E0,I0,F0,q0);
  subplot(1,4,1);
  plot(x,Ut,'o-',x,v); legend("solution MEF","solution exacte"); grid(); title("deplacement pour un maillage variable (m)");
  % de meme pour le maillage regulier
  x=maillage_regulier(nh,l);
  q=x*q0/l;
  [k,f]=assemblage(nh,l,E,I,q,0)
  u=dep(nh,l)
  Q=efforts(nh,l,F0)
  u=resolutionU(u,k,Q+f)
  % pour dessiner le graphique on prend uniquement les valeurs de u relatives a des translations
  j=1;
  for i=1:2*nh+2
    if mod(i,2)==1
      Ut(j)=u(i);
      j=j+1;
    endif
  endfor
  v=vAnalytique(x,l,E0,I0,F0,q0);
  subplot(1,4,2);
  plot(x,Ut,'o-',x,v); legend("solution MEF","solution exacte"); grid(); title("deplacement pour un maillage regulier (m)");
  Qele1=Qele(1,u,x,E,I,q)
  printf("Ty(0)=%f N | Mfz(0)=%f Nm | Ty(he)=%f N | Mfz(he)=%f Nm \n",-Qele1(1),-Qele1(2),Qele1(3),Qele1(4));
 
  % valeurs exactes de Ty Mfz
  Ty=F0+q0*l/2*(1-x.^2/l^2)                    
  Mfz=F0*(l-x)+q0*l^2/6*(2-3*x/l+x.^3/l^3) 
  
  % valeurs par la MEF
  for e=1:nh                                            % il y a des valeurs qui se repetent
    Ty_mef(e)=-Qele(e,u,x,E,I,q)(1);                    % pour des elements consecutifs
  endfor                                                % => on prend une valeur par element
  Ty_mef(nh+1)=Qele(nh,u,x,E,I,q)(3); 
  Ty_mef
  
  for e=1:nh
    Mfz_mef(e)=-Qele(e,u,x,E,I,q)(2);
  endfor
  Mfz_mef(nh+1)=Qele(nh,u,x,E,I,q)(4);
  Mfz_mef
  
  subplot(1,4,3);
  plot(x,Ty_mef,'o-',x,Ty); legend("solution MEF","solution exacte"); grid(); title("Ty (N)");
  subplot(1,4,4);
  plot(x,Mfz_mef,'o-',x,Mfz); legend("solution MEF","solution exacte"); grid(); title("Mfz (Nm)");
  
endfunction
main()

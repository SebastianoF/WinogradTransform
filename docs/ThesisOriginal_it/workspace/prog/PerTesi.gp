



\\\/////////////////////////////
\\\/////////////////////////////
\\   IMPLEMENTAZIONI ESSENZIALI
\\\/////////////////////////////
\\\/////////////////////////////


\\SHIFTERS

\\Shifter nell'algebra dei vettori del vettore a
\\ input: vettore a
\\ output: shift a destra del vettore a
{
vectShifter(a) = local(r,b,ris);
  r = length(a);
  b = vecextract(a,concat("1..",Str(n-1)));
  ris = concat(a[n],b);
ris;

}

\\Shifter nell'algebra dei vettori del vettore a
\\ input: vettore a
\\ output: shift a destra del vettore a
{
vecShifterK(a,k) = local(ans);
  ans=a;
  for(i = 1, k,
    ans = vecShifter(ans)
  );
ris;
}


\\shifter dell'algebra delle matrici
\\input: intero positivo
\\output: shifter r-dimensionale dell'algebra delle matrici
{
matShifter(r) = local(ans);
  ans = matrix(r,r,i,j,
          if(i == j-1,1,0);
        );
ans;
}


\\ TOOLS


{
\\elimina le ripetizioni nel vettore v
elimina(v)=local();
  eval(Set(v));
}


{
\\dati due vettori a,b restituisce il vettore (senza ripetizioni) formato dagli elementi di a che non sono in b
complemento(a,b)=local();
  eval(setminus(Set(a),Set(b)));
}


\\ Input: n, gruppo H contenuto in (Z_n)*
\\ Output: insieme delle orbite di Z_n sotto l'azione del gruppo H contenuto in (Z_n)*
\\ Usa: elimina, complemento
{
classi(n,gruppo)=local(nn,X,G,ris,gg,ele,orbita);
  X=vector(n,k,k-1);
  gg=length(gruppo);
  if(gg==0,return(0),G=gruppo);
  ris=[];
  nn=n;
  while(nn!=0,ele=Mod(X[1],n);
  orbita=[];
  for(ii=1,gg,orbita=concat(orbita,G[ii]*ele));
  orbita=elimina(lift(orbita));
  ris=concat(ris,[orbita]);
  X=complemento(X,orbita);
  nn=length(X));
ris;
}



{
\\restituisce il gruppo generato da q in (Z_n)* (q deve essere coprimo con n)
generato(n,q)=local(g,ris,e);
  if(gcd(q,n)!=1,return(0));
   g=Mod(q,n);ris=[];
   for(k=1,n-1,
      e=lift(g^k);
      ris=concat(ris,e);
      if(e==1,return(ris))
   );
}


\\inseieme delle etichette L
\\ Input: interi r, q = p^n
\\ Output: vettore L dei gradi dei polinomi irriducibili della fattorizzazione x^r - 1 in GF(q)[x]
\\ usa: classi, generato
{
vettoreL(n,q)=local(z);
  z=classi(n,generato(n,q));
  s=length(z);
  v=vector(s,j,length(z[j]));
vecsort(v);
}


{
\\trova i gradi dei polinomi irriducibili
\\nel cui prodotto si fattorizza x^n-1 in GF(q)[x]
gradi(n,q)=local(z);
  z=classi(n,generato(n,q));
  s=length(z);
  v=vector(s,j,length(z[j]));
vecsort(v);
}

\\ TRASFORMAZIONI PSI

\\ Input: vettore a
\\ Output: matrice circolante m definita dal vettore a
{
psi1(a) = local(r,s,m);
  r = length(a);
  s = matShifter(r);
  m = matrix(r,r);
  for(i=0, r-1,
    m = m + a[i+1]*s^i
  );
m;
}

\\ Input: matrice circolante definita dal vettore a
\\ Output: vettore a
{
psi1_inv(m) = local(a);
  a = m[1,];
a;
}


\\ Input: polinomio a = a(x) di R_{r,q}
\\ Output: vettore corrispondente di dimensione r
{
psi2(a,r) = local(v);
  a = lift(lift(a));
  v = concat(Vec(a) , vector(r - poldegree(a)));
v;
}


\\ Input: vettore v di dimensione r, intero r
\\ Output: polinomio in R_{r,q} definito da v
{
psi2_inv(v,r) = local(a);
  a = Polrev(v)*Mod(1, x^r - 1);
a;
}


\\ Input: u elemento dell'algebra FC_{r} con g
\\        elemento del gruppo ciclico,
\\        r intero postivo
\\ Output: matrice circolante corrispondente
{
psi3(u,r) = local(v, m);
  v = concat(Vec(u) , vector(r - poldegree(u)));
  m = psi1(v);
m;
}


\\ Input:  matrice circolante m
\\ Output: elemento dell'algebra FC_{r} corrisp.
{
psi3_inv(m) = local(r, a);
  r = matsize(m)[1];
  for(i = 0, r-1,
    a = m[1,i+1]*g^(i)
  );
a;
}



\\psi5
\\ Input: a vettore di polinomi in di Q_{r,q}
\\ Output: elemento dell'algebra FC_{r} corrisp
{
psi5(u,r,q) = local(vett, fc);
  fc = 'da fare';
fc;
}


/////////////////////////////////////
////////// Winograd

{
\\resto della divisione di a(x) per b(x) in Zp[x]
restop(a,b,p)=lift(divrem(Mod(1,p)*a,b)[2])
}


{
\\tcrp : riceve un array di polinomi, un array dei relativi moduli e il modulo p e restituisce il polinomio soluzione
tcrp(ap,am,p)=local(n,u,v);
   n=length(ap);
   v=vector(n,k,Mod(Mod(1,p)*ap[k],am[k]));
   u=chinese(v);
return(lift(lift(u)))
}



{
\\quoziente della divisione di a(x) per b(x) in Zp[x]
quozp(a,b,p)=lift(divrem(Mod(1,p)*a,b)[1])
}


\\isomorfismo tra Rp,n e il relativo prodotto di campi.
\\riceve un polinomio in Rp,n (ovvero un polinomio di grado < n con coeff. interi, che sono visti modulo p)
\\restituisce un array di s polinomi, dove s è il numero dei poli irriducibili nei quali x^n-1 si fattorizza modulo p
{
teta(a,n,p)=local(f,s);
  f=lift(factormod(x^n-1,p));
  s=matsize(f)[1];
return(vector(s,i,restop(a,f[i,1],p)))
}

{
\\inverso di teta
gamm(a,n,p)=local();
   tcrp(a,circolmax(n,p),p);
}



{
\\restituisce la matrice T che rappresenta tetaQ
del(n)=local(ff,poli,s,T,resti);
ff=factor(x^n-1);
poli=ff[,1]~;
s=length(poli);
T=matrix(n,n);resti=vector(n);
for(j=1,n,resti[j]=[];for(i=1,s,resti[j]=concat(resti[j],pol2vetc(divrem(x^(j-1),poli[i])[2],poldegree(poli[i])))));
for(i=1,n,T[,i]=mattranspose(resti[i]));
T;
}

{
\\applica del(n) come pura trasformazione lineare
dellin(a)=local(n);
n=length(a);
mattranspose(del(n)*mattranspose(a));
}

{
\\applica del(n) come biiezione tra Q[x]/(x^1-1) e il prodotto dei Q[x]/(polciclotomico fattore di x^n-1)
\\è la formattazione corretta di dellin(a,n)
\\dà lo stesso risultato di tetaQ(a,n)
delpol(a,n)=local(ap,b,ff,poli,s,v,d);
ap=pol2vetc(a,n);
b=dellin(ap);
ff=factor(x^n-1);
poli=ff[,1]~;
s=length(poli);
v=vector(s);ini=0;
for(k=1,s,d=poldegree(poli[k]);v[k]=Polrev(vector(d,j,b[ini+j]));ini+=d);
v;
}

{
\\inversa di del(n)
gamlin(a)=local(n,mm);
n=length(a);
mm=del(n)^-1;
mattranspose(mm*mattranspose(a));
}

{
\\dà lo stesso risultato di del(n)^-1, ma costruisce la matrice senza bisogno di conoscere del(n)
gam(n)=local(ff,poli,s,ris,ma,d,ini);
ff=factor(x^n-1);
poli=ff[,1]~;
s=length(poli);
ris=matrix(n,n);ini=1;
ma=matid(s);
for(k=1,s,d=poldegree(poli[k]);ris[ini,]=pol2vetc(tcrQ(ma[k,],poli),n);for(j=ini+1,ini+d-1,ris[j,]=destra(ris[j-1,]));ini+=d);
mattranspose(ris);
}

{
\\restituisce la matrice T che rappresenta teta
delnp(n,p)=local(ff,poli,s,T,resti);
if(gcd(n,p)!=1,return(0));
ff=lift(factormod(x^n-1,p));
poli=ff[,1]~;
s=length(poli);
T=matrix(n,n);resti=vector(n);
for(j=1,n,resti[j]=[];for(i=1,s,resti[j]=concat(resti[j],pol2vetc(restop(x^(j-1),poli[i],p),poldegree(poli[i])))));
for(i=1,n,T[,i]=mattranspose(resti[i]));
T;
}

{
\\dà lo stesso risultato di delnp(n,p)^-1, ma costruisce la matrice senza bisogno di conoscere delnp(n,p)
gamnp(n,p)=local(ff,poli,s,ris,ma,d,ini);
ff=lift(factormod(x^n-1,p));
poli=ff[,1]~;
s=length(poli);
ris=matrix(n,n);ini=1;
ma=matid(s);
for(k=1,s,d=poldegree(poli[k]);ris[ini,]=pol2vetc(tcrp(ma[k,],poli,p),n);for(j=ini+1,ini+d-1,ris[j,]=destra(ris[j-1,]));ini+=d);
mattranspose(ris);
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\/////////////////////////////////////////////////////////////////
\\ idempotenti


{
\\idempotenti minimali
circolidem(n,p)=local(c,id,s,ris);
c=circolmax(n,p);
s=length(c);
ris=vector(s);
id=matid(s);
for(j=1,s,ris[j]=tcrp(id[j,],c,p));
ris;
}

{
\\decomposizione nella somma
decsomma(a,n,p)=local(ee,s);
ee=circolidem(n,p);s=length(ee);
vector(s,k,moltnp(a,ee[k],n,p));
}

{
\\decomposizione nel prodotto
decprod(a,n,p)=local(ee,s);
ee=circolidem(n,p);s=length(ee);
vector(s,k,moltnp(1,1+moltnp(ee[k],a-1,n,p),n,p));
}

{
\\decomposizione di una matrice circolante nella somma (modulo p)
decmatsomma(ma,p)=local(ee,s,a,n,vv);
a=Polrev(ma[1,]);n=matsize(ma)[1];
ee=circolidem(n,p);s=length(ee);
vv=vector(s,k,moltnp(a,ee[k],n,p));
vector(s,k,pol2mat(vv[k],n));
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\/////////////////////////////////////////////////////////////////

{
\\decomposizione di una matrice circolante nel prodotto (modulo p)
decmatprod(ma,p)=local(ee,s,a,n,vv);
a=Polrev(ma[1,]);n=matsize(ma)[1];
ee=circolidem(n,p);s=length(ee);
vv=vector(s,k,moltnp(1,1+moltnp(ee[k],a-1,n,p),n,p));
vector(s,k,pol2mat(vv[k],n));
}




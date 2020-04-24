
/*

LABORATORIO DI CALCOLO SIMBOLICO
ANNO ACCADEMICO 2010/2011

SEBASTIANO FERRARIS

\r Ferraris\lezione_12

*/

\\____DODICESIMA LEZIONE_________________________________________________
\\____________________________________________________________________
\\____________________________________________________________________

\\____________________________________________________________________
\\Programmi dati in classe LEZ 12

{\\restituisce la lista dei gap tra i primi compresi tra a e b
gaps(a,b)=local(lis,gg,nn);
lis=[];forprime(n=a,b,lis=concat(lis,n));
nn=length(lis);gg=[];
for(i=1,nn-1,gg=concat(gg,lis[i+1]-lis[i]));
gg;
}

{\\restituisce la lista delle posizioni di x nella lista lista
posizioni(lista,x)=local(nn,ris);
nn=length(lista);ris=[];
for(kk=1,nn,if(lista[kk]==x,ris=concat(ris,kk)));
ris;
}

{\\restituisce il numero dei gap di lunghezza m tra a e b
gaps1(a,b,m)=local(lis,gg,nn,ss,cont,ris);
lis=gaps(a,b);
length(posizioni(lis,m));
}

{\\restituisce il primo successivo a n (vedi nextprime)
primodopo(n)=local(aa);
aa=n+1;
while(!isprime(aa),aa++);
aa;
}

{\\riceve una lista di interi e restituisce la stessa con l'aggiunta del minimo primo che divide il prodotto della lista + 1
euc1(lista)=local(qr,h);
qr=prod(i=1,length(lista),lista[i]);
h=factor(qr++);
concat(lista,h[1,1]);
}

{\\usa euc1 n volte e scrive il risultato nel file nomefile ES: useuc1([2],16,"cerruti/eu16")
useuc1(lista,n,nomefile)=local(a);
a=lista;
for(i=1,n,a=euc1(a));
write(nomefile,a)}

{\\schema dell'algoritmo potenza veloce; operid è l'elemento neutro, oper è l'operazione
pot(a,m)=local(u,t,v);
u = m; t = operid; v = a;
while(u != 0, if(u%2 == 1, t = oper(t,v)); v = oper(v,v);
u = floor(u/2));
t;
}

{\\shift a destra di un posto
destra(a)=local(n,b,ris);
n=length(a);b=vecextract(a,concat("1..",Str(n-1)));ris=concat(a[n],b);
ris;
}

{\\shift a destra di k posti
destrak(a,k)=local(ris);
ris=a;for(i=1,k,ris=destra(ris));
ris;
}

{\\shift a sinistra di un posto
sinistra(a)=local();
destrak(a,length(a)-1);
}

{\\shift a sinistra di k posti
sinistrak(a,k)=local(ris);
ris=a;for(i=1,k,ris=sinistra(ris));
ris;
}

{\\test di Fermat
fermatest(b,n)=local(g);
g=Mod(b,n);
if(lift(g^(n-1))==1,1,0);
}

{\\lista degli pseudoprimi di Fermat
fermatlista(b,n)=local(g);
g=[];
forstep(k=9,n,2,if(fermatest(b,k)==1 && isprime(k)==0,g=concat(g,k)));
g;
}

{\\test di Eulero
eulertest(b,n)=local(g,h);
g=Mod(b,n);h=Mod(kronecker(b,n),n);
if(g^((n-1)/2)==h,1,0);
}

{\\lista degli pseudoprimi di Eulero
eulerlista(b,n)=local(g);
g=[];
forstep(k=9,n,2,if(eulertest(b,k)==1 && gcd(b,k)==1 && isprime(k)==0,g=concat(g,k)));
g;
}

{\\lista dei generatori di n
generatori(n)=local(b,e,ris,phi);
if(n<2,return(0));if(n==2,return([1]));if(n==4,return([3]));ris=[];phi=eulerphi(n);
e=ispower(n);
if(e>0,b=floor(n^(1/e));if(b==2 || !isprime(b), return([]));for(i=2,n-2,if(gcd(i,n)==1 && znorder(Mod(i,n))==phi,
ris=concat(ris,i)));return(ris));
if(!isprime(n),return([]),for(i=2,n-2,if(znorder(Mod(i,n))==phi,ris=concat(ris,i))));
ris;
}

{\\riceve un array intmod e restituisce la soluzione
cinese(lista)=local(n,u,v);
n=length(lista);u=chinese(lista[1],lista[2]);
for(i=3,n,v=chinese(u,lista[i]);u=v);
u;
}

{\\risolve il sistema {x=2 (mod 3), x=3 (mod 5), ..., x=j (mod p_j), ... } per 1<j<=n
\\ sequenza A053664 in OEIS
gioco(n)=local(vv,ris);
ris=[];
for(k=2,n,vv=vector(k,x,Mod(x,prime(x)));ris=lift(concat(ris,cinese(vv))));
ris;
}

{\\restituisce 1 se n è il periodo di elem (con operazione oper e neutro operid, vedi pot)
periodo(elem,n)=local(ff,s);
if(pot(elem,n)!=operid,return(0));
ff=factor(n);
s=matsize(ff)[1];
for(i=1,s,if(pot(elem,n/ff[i,1])==operid,return(0)));
1;
}

{\\trova il periodo di elem (con operazione oper e neutro operid, vedi pot) se ordine è un suo multiplo
cercaperiodo(elem,ordine)=local();
fordiv(ordine,n,if(pot(elem,n)==operid,return(n)));
0;
}

{
diric1(n,lista)=local(g,lp,pr,qr,pol,t,h);
g=factor(n);
lp=g[,1]~;pr=prod(i=1,length(lp),lp[i]);
qr=prod(i=1,length(lista),lista[i]);pol=polcyclo(n);
t=1;while(abs(substpol(pol,x,t*pr*qr))<2,t=t+1);
h=factor(abs(substpol(pol,x,t*pr*qr)));
return(concat(lista,h[,1]~))}


{
diric1min(n,lista)=local(g,lp,pr,qr,pol,t,h);
g=factor(n);
lp=g[,1]~;pr=prod(i=1,length(lp),lp[i]);
qr=prod(i=1,length(lista),lista[i]);pol=polcyclo(n);
t=1;while(abs(substpol(pol,x,t*pr*qr))<2,t=t+1);
h=factor(abs(substpol(pol,x,t*pr*qr)));
return(concat(lista,h[1,1]))}


{
mersenne(p)=local(n,a);
n=2^p-1;a=Mod(4,n);
for(i=1,p-2,a=a^2-2);
if(lift(a)==0,return("OK"),return("NO"))}

{
lucas(n,lim)=local(qf,zz,le,ris);
qf=factor(n-1)[,1];zz=vector(lim,k,prime(k));
le=length(qf);
for(i=1,le,qq=qf[i];ris=lucascond(n,qq,zz);
if(ris==0,print(qq);return(0));
if(ris==-1,return(-1)));
return(1)}


{
lucascond(n,qq,zz)=local(le,a,b,c);le=length(zz);
for(i=1,le,a=Mod(zz[i],n);b=lift(a^(n-1));c=lift(a^((n-1)/qq));
if(b!=1,return(0));
if(c!=1,return(1)));
return(-1)}

{
pock(n,F,base)=local(qF,ris,le);
qF=factor(F)[,1];
le=length(qF);
for(i=1,le,qq=qF[i];
ris=pockcond(n,qq,base);
if(ris==0,return(0));
if(ris==-1,return(-1)));
return(1)}


{
pockcond(n,qq,base)=local(a,b,c,d);
a=Mod(base,n);
b=lift(a^(n-1));
c=lift(a^((n-1)/qq));
d=gcd(c-1,n);
if(b!=1,return(0));
if(d==1,return(1),return(-1));
}

{
cercapock(F,base,ini,fine)=local(ris);
ris=[];
for(i=ini,fine,n=i*F+1;if(pock(n,F,base)==1,ris=concat(ris,n)));
return(ris)}


{
proth(m,k,base)=local(n,a,c);
if(2^k<m,return(-2));
n=m*2^k+1;
if(gcd(n,base)!=1,return(0));
a=Mod(base,n);
c=lift(a^((n-1)/2));
if(c==n-1,return(1),-1)}


{
cercaproth(k,base)=local(ris);
ris=[];
for(m=1,2^k-1,if(proth(m,k,base)==1,ris=concat(ris,m*2^k+1)));
return(ris)}

{
akstest(a,r,n)=local();
if(gcd(n,r)!=1,return(0));
if(gcd(n,a)!=1,return(0));
xx=Mod(x,x^r-1);
bb=Mod(a,n);
return((xx+bb)^n==xx^n+bb);
}

{
akslista(a,r,n)=local(g);
g=[];
forstep(k=9,n,2,if(akstest(a,r,k) && !isprime(k),g=concat(g,k)));
g;
}

{\\composto di permutazioni
\\(2,3,1) significa 1->2, 2->3, 3->1 etc.
permcomp(a,b)=local(n,c);
n=length(a);
c=vector(n);
for(i=1,n,c[i]=a[b[i]]);
c;
}

{\\permutazione identica
permid(n)=local();
vector(n,j,j);
}

{\\matrice di permutazione associata alla permutazione a
permat(a)=local(n,c,w,v);
n=length(a);
c=matrix(n,n);
w=[i,j];v=[[1,n],[1,n]];forvec(w=v,if(a[w[2]]==w[1],c[w[1],w[2]]=1));
c;
}

{\\inversa di permat
permatinv(a)=local(ris);
ris=vector(matsize(a)[1]);
for(i=1,n,ris[i]=posizioni(a[,i],1)[1]);
ris;
}

{\\matrice di shift a destra
matshiftd(n)=local();
permat(destra(permid(n)));
}

{\\matrice circolante destra con prima riga a
circol(a)=local(n,s,ris);
n=length(a);
s=matshiftd(n);
ris=matrix(n,n);
for(i=0,n-1,ris=ris+a[i+1]*s^i);
ris;
}

{\\inverte l'ordine in un vettore
contr(v)=local(ris,n);
n=length(v);ris=vector(n);forstep(x=n,1,-1,ris[n-x+1]=v[x]);
ris;
}

\\vettore -> polinomio
\\usare Polrev(vettore, {var})

{\\polinomio -> vettore
\\inverso di Polrev
pol2vet(p)=local(ris);
ris=Vec(p);contr(ris);
}

{
\\polinomio -> vettore di lunghezza data
pol2vetc(p,n)=local(ris);
ris=Vec(p);m=length(ris);if(n<m,return(0),return(contr(concat(vector(n-m),ris))));
}

{\\riceve due polinomi a,b in x, di grado <= n-1 e calcola a*b mod x^n-1 , dove n è un intero positivo
moltpoln(a,b,n)=local(y,z);
y=Mod(a,x^n-1);z=Mod(b,x^n-1);
lift(y*z);
}

{\\restituisce l'elenco dei cicli disgiunti della permutazione g
cicli(g)=local(nn,ciclo,ris,s,ini);
nn=length(g);
s=vector(nn,n,n);ris=[];
while(length(s)>0,ini=s[1];x=ini;ciclo=[x];while(g[x]!=ini,ciclo=concat(ciclo,g[x]);x=g[x]);s=complemento(s,ciclo);ris=concat(ris,[ciclo]));
ris;
}

{\\calcola il periodo di g
permper(g)=local(z,v);
z=cicli(g);v=vector(length(z),i,length(z[i]));
lcm(v);
}

{\\calcola il periodo di g (permutazione di 1..n), con cercaperiodo; si deve porre esternamente oper(a,b)=permcomp(a,b) e operid=permid(n)
permper2(g)=local(z,v);
cercaperiodo(g,lcm(operid));
}

{\\matrice circolante -> polinomio che la rappresenta
mat2pol(a)=local();
Polrev(a[1,]);
}

{\\inversa di mat2pol
pol2mat(b,n)=local(c);
c=pol2vetc(b,n);
circol(c);
}

{\\radice primitiva n-esima modulo p primo; richiede che n divida p-1
radnp(n,p)=local(z);
if(gcd(n,p)!=1||!isprime(p),return(0));
if((p-1)%n!=0,return(0));
z=znprimroot(p);
z^((p-1)/n);
}

{\\matrice di Fourier nxn con n che divide p-1, p primo
fournp(n,p)=local(z);
z=radnp(n,p);
matrix(n,n,i,j,z^((i-1)*(j-1)));
}

{\\inversa di Fourier; fournp(n,p)*Mod(circol(vettore),p)*fournpinv(n,p) è diagonale
fournpinv(n,p)=local(z,m,s);
z=radnp(n,p);
m=matrix(n,n,i,j,z^(-(i-1)*(j-1)));
s=Mod(n,p);
s^(-1)*m;
}

{\\restituisce W_n(a,b,h,k), utilizzando i polinomi
ric(a,b,h,k,n)=local(xx,zz);
xx=Mod(X,X^2-h*X+k);
xx=xx^n;
zz=Vec(lift(xx));
if(#zz==1,zz=[0,zz[1]]);
return(b*zz[1]+a*zz[2])}

{\\dà ric(a,b,h,k,n) modulo m, utilizzando la matrice compagna
ricmod(a,b,h,k,n,m)=local(z,y);
z=Mod([0,1;-k,h],m);
y=z^n;
return(lift(a*y[1,1]+b*y[1,2]))}

{\\funzione psidelta, dove delta=h^2-4k
psidelta(h,k,n)=local(delta,ff);
delta=h^2-4*k;if(gcd(n,delta*k)>1,return(0));
ff=factor(n);
prod(k=1,matsize(ff)[1],ff[k,1]^(ff[k,2]-1)*(ff[k,1]-kronecker(delta,ff[k,1])));
}

{\\prodotto nel gruppo G_{n,x^2-hx+k}
gprod(h,k,n,u,v)=local(u1,v1,r,ve);
u1=Mod(Mod(1,n)*u,x^2-h*x+k);v1=Mod(Mod(1,n)*v,x^2-h*x+k);
r=lift(lift(u1*v1));
ve=Vec(r);
if(length(ve)==1,return(1));
d=gcd(ve[1],ve[2]);
(ve[1]/d)*x+ve[2]/d;
}

{\\elimina le ripetizioni nel vettore v
elimina(v)=local();
eval(Set(v));
}

{\\dati due vettori a,b restituisce il vettore (senza ripetizioni) formato dagli elementi di a che non sono in b
complemento(a,b)=local();
eval(setminus(Set(a),Set(b)));
}

{\\restituisce l'insieme dei coprimi con n
coprimi(n)=local();
ris=[];
for(i=1,n-1,if(gcd(i,n)==1,ris=concat(ris,i)));
ris;
}

{\\restituisce il gruppo generato da q in (Z_n)* (q deve essere coprimo con n)
generato(n,q)=local(g,ris,e);
if(gcd(q,n)!=1,return(0));
g=Mod(q,n);ris=[];
for(k=1,n-1,e=lift(g^k);ris=concat(ris,e);if(e==1,return(ris)));
}

{\\calcola l'insieme delle orbite di Z_n sotto l'azione del gruppo H contenuto in (Z_n)*
classi(n,gruppo)=local(nn,X,G,ris,gg,ele,orbita);
X=vector(n,k,k-1);gg=length(gruppo);if(gg==0,return(0),G=gruppo);ris=[];nn=n;
while(nn!=0,ele=Mod(X[1],n);orbita=[];
for(ii=1,gg,orbita=concat(orbita,G[ii]*ele));
orbita=elimina(lift(orbita));ris=concat(ris,[orbita]);X=complemento(X,orbita);nn=length(X));
ris;
}

{\\risolve, se possibile, x^2 = a (mod p)
qsolve(a,p)=local(unmezzo,h,k);
if(gcd(a,p)!=1,return(0));
if(!isprime(p) || kronecker(a,p)!=1,return(0));
if(p==2,return(1));
unmezzo=(p+1)/2;
h=1;k=a;
while(kronecker(h^2-4*k,p)!=-1,h=h+1);
ris=ricmod(2,h,h,k,(p+1)/2,p);
ris=(unmezzo*ris)%p;
ris=min(ris,p-ris);
return(ris);}

{matinvQ(a)=local(n);
n=matsize(a)[1];
if(gcd(Polrev(a[1,]),x^n-1)==1,1,0);
}

{
apparizione(h,k,n)=local(delta,u,t,v);
delta=h^2-4*k;if(gcd(n,delta*k)>1,return(0));
fordiv(psidelta(h,k,n),m,	
u = m; t = 1; v = x;
while(u != 0, if(u%2 == 1, t = gprod(h,k,n,t,v)); v = gprod(h,k,n,v,v);
u = floor(u/2));
if(t==1,return(m)));
}

{
gradi(n,q)=local(z);
z=classi(n,generato(n,q));
s=length(z);
v=vector(s,j,length(z[j]));
vecsort(v);
}

{
circolgrorder(n,q)=local(v);
v=gradi(n,q);
prod(i=1,length(v),q^(v[i])-1);
}

{
factoringfq(a,p,e)=local(f);
f=lift(ffinit(p,e,t));
lift(lift(factorff(a,p,f)));
}

{
circolmax(n,p)=local();
(lift(factormod(x^n-1,p))[,1])~
}

{
circolmin(n,p)=local(g);
g=factormod(x^n-1,p);
vector(matsize(g)[1],k,lift((x^n-1)/g[k,1]));
}

{
\\riceve tre polinomi a,b,p in una stessa variabile e calcola a*b mod (p,n) , dove n è un intero positivo
moltp(a,b,p,n)=local(y,z);
y=Mod(Mod(1,n)*a,p);z=Mod(Mod(1,n)*b,p);
return(lift(lift(y*z)))}

{
\\riceve due polinomi in una stessa variabile e due interi e=esponente, n=modulo
\\n deve essere positivo, e può essere negativo, se a è coprimo con p
\\calcola a^e mod (p,n)
potp(a,e,p,n)=local(y);
y=Mod(Mod(1,n)*a,p);
return(lift(lift(y^e)))}

{
\\calcola il mcd tra due polinomi a, b in Zp[x]
mcdp(a,b,p)=local(ap,bp);
ap=Mod(1,p)*a;bp=Mod(1,p)*b;
return(lift(gcd(ap,bp)))}

{
\\moltiplicazione di polinomi nell'algebra Rp,n = Zp[x]/(x^n-1)
moltnp(a,b,n,p)=moltp(a,b,x^n-1,p)}


{
\\potenza di un polinomio nell'algebra Ap,n
potnp(a,e,n,p)=potp(a,e,x^n-1,p)}


{
\\tcrp : riceve un array di polinomi, un array dei relativi moduli e il modulo p e restituisce il polinomio soluzione
tcrp(ap,am,p)=local(n,u,v);
n=length(ap);v=vector(n,k,Mod(Mod(1,p)*ap[k],am[k]));u=chinese(v);
return(lift(lift(u)))}

{
\\resto della divisione di a(x) per b(x) in Zp[x]
restop(a,b,p)=lift(divrem(Mod(1,p)*a,b)[2])}


{
\\quoziente della divisione di a(x) per b(x) in Zp[x]
quozp(a,b,p)=lift(divrem(Mod(1,p)*a,b)[1])}

{
\\isomorfismo tra Rp,n e il relativo prodotto di campi.
\\riceve un polinomio in Ap,n (ovvero un polinomio di grado < n con coeff. interi, che sono visti modulo p)
\\restituisce un array ci s polinomi, dove s è il numero dei poli irriducibili nei quali x^n-1 si fattorizza modulo p
teta(a,n,p)=local(f,s);
f=lift(factormod(x^n-1,p));s=matsize(f)[1];return(vector(s,i,restop(a,f[i,1],p)))}

{
\\inverso di teta
gamm(a,n,p)=local();
tcrp(a,circolmax(n,p),p);
}

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

{
\\decomposizione di una matrice circolante nel prodotto (modulo p)
decmatprod(ma,p)=local(ee,s,a,n,vv);
a=Polrev(ma[1,]);n=matsize(ma)[1];
ee=circolidem(n,p);s=length(ee);
vv=vector(s,k,moltnp(1,1+moltnp(ee[k],a-1,n,p),n,p));
vector(s,k,pol2mat(vv[k],n));
}

{
\\tcrQ : riceve un array di polinomi, un array dei relativi moduli, e restituisce il polinomio soluzione nei razionali
tcrQ(ap,am)=local(n,u,v);
n=length(ap);v=vector(n,k,Mod(ap[k],am[k]));u=chinese(v);
return(lift(u));
}

{
\\idempotenti minimali razionali
circolidemQ(n)=local(f,c,id,s,ris);
f=factor(x^n-1);c=f[,1]~;
s=length(c);
ris=vector(s);
id=matid(s);
for(j=1,s,ris[j]=tcrQ(id[j,],c));
ris;
}

{
\\moltiplica i polinomi modulo x^n-1
moltnQ(a,b,n)=local();
divrem(a*b,x^n-1)[2];
}

{
\\decomposizione di una matrice circolante nella somma (su Q)
decmatsommaQ(ma)=local(ee,s,a,n,vv);
a=Polrev(ma[1,]);n=matsize(ma)[1];
ee=circolidemQ(n);s=length(ee);
vv=vector(s,k,moltnQ(a,ee[k],n));
vector(s,k,pol2mat(vv[k],n));
}

{
\\decomposizione di una matrice circolante nel prodotto (su Q)
decmatprod(ma)=local(ee,s,a,n,vv);
a=Polrev(ma[1,]);n=matsize(ma)[1];
ee=circolidemQ(n);s=length(ee);
vv=vector(s,k,1+moltnQ(ee[k],a-1,n));
vector(s,k,pol2mat(vv[k],n));
}

{
\\isomorfismo tra Q[x]/(x^n-1) e il relativo prodotto di campi.
tetaQ(a,n)=local(f,s);
f=factor(x^n-1);s=matsize(f)[1];
vector(s,i,divrem(a,f[i,1])[2]);
}

{
\\inverso di tetaQ
gammaQ(a,n)=local(f,c);
f=factor(x^n-1);c=f[,1]~;
tcrQ(a,c);
}

{
\\vettore casuale di lunghezza n con componenti comprese tra 0 e m-1
randvett(n,m)=vector(n,k,random(m));
}

{
\\probabilità che una circolante nxn abbia il determinante divisibile per il primo p, n generico
npprobgen(n,p)=local();
if(n%p==0,return(0));
1.-circolgrorder(n,p)/(p^n);
}

{
\\probabilità che una circolante nxn abbia il determinante divisibile per il primo p, posto che n divida p-1
npprob(n,p)=local();
if((p-1)%n != 0,return(0));
1.-(1-1/p)^n;
}

{
\\probabilità ricavata sperimentalmente
npprobex(n,p,casi)=local();
favorevoli=0;
for(i=1,casi,if(matdet(circol(randvett(n,p)))%p==0,favorevoli++));
1.*favorevoli/casi;
}

{
\\log in base 2
log2(a)=log(a)/log(2);
}

{
\\riceve n
\\restituisce un fattore di n
\\o un r coprimo con r tale che il periodo moltiplicativo di n mod r sia > [4log(n)^2]
akspasso2(n)=local(r,L,v);
r=2;L=floor(4*log2(n)^2);
while(r<n,d=gcd(r,n);if(d>1,return(d));v=znorder(Mod(n,r));
if(v>L,break());
r=r+1);
[r,v,L,isprime(r),n,floor(16*log(n)^5)]
}

{\\AKS
aks(n)= local(w,r,L,l,d,v,xx,bb);
if(n==2,return("primo"));
\\termina se n è pari diverso da 2
if(n%2==0,print("n è pari");return(0));
\\verifica se è una potenza perfetta
w=ispower(n);
if(w>0,print("n è una potenza perfetta");return(w));
\\trova r
r=2;L=floor(4*(log2(n))^2);
while(r<n,d=gcd(r,n);if(d!=1,print("n è composto");return(d));v=znorder(Mod(n,r));
if(v>L,break());
r=r+1);
if(r==n, return("primo"));
\\loop finale
l=2*(floor(sqrt(r)))*floor(log2(n))+1;
xx=Mod(x,x^r-1);
for(b=1,l,
bb=Mod(b,n);
if((xx+bb)^n!=xx^n+bb,print("n è composto");return(bb)));
return(1);
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

{
provaecm(nn,ee)=local(xx,yy,aa,bb,elli);
xx=random(nn);yy=random(nn);aa=random(nn);
bb=yy^2-xx^3-aa*xx;
elli=ellinit([0,0,0,Mod(aa,nn),Mod(bb,nn)]);
ellpow(elli,[xx,yy],ee);
}


{
hadamard(n)=local(hh,mh,m1,m2);
hh=[1,1;1,-1];
if(n==1,return(hh));mh=hh;
for(j=1,n-1,m1=concat(mh,mh);m2=concat(mh,-mh);mh=concat(mattranspose(m1),mattranspose(m2)));
mh;
}


{
	ellido(nn)=local(xx,yy,aa,bb);
	xx=random(nn);yy=random(nn);aa=random(nn);
	bb=yy^2-xx^3-aa*xx;
	[ellinit([0,0,0,Mod(aa,nn),Mod(bb,nn)]),[xx,yy]];
}


{
ellisum(a,P,Q)=local(lam, mu, x3, x1, y1, x2, y2);
 if(P == [0], return(Q)); if(Q == [0], return(P)); x1 = P[1]; x2 = Q[1]; y1 = P[2]; 
 y2 = Q[2];if(x1==x2&&y1==-y2,return([0]));
 if(P == Q, lam = (3*x1^2 + 2*a[2]*x1 + a[4] - a[1]* y1)/(2* y1 + a[1]* x1 + a[3]), lam = (y2 - y1)/(x2 - x1)); 
 mu = y1 - lam *x1; x3 = lam^2 + a[1]* lam - a[2] - x1 - x2;
 [x3, -(lam + a[1])* x3 - mu - a[3]];
}


{
ellipot(a,P,m)=local(u,t,v);
u = m; t = [0]; v = P;
while(u != 0, if(u%2 == 1, t = ellisum(a,t,v)); v = ellisum(a,v,v);
u = floor(u/2));
t;
}


{
ellisumn(z,P,Q,n)=local(lam, mu, x3, x1, y1, x2, y2,den,fat,a);
 if(P == [0], return(Q)); if(Q == [0], return(P)); 
 a=z*Mod(1,n);
 x1 = Mod(P[1],n); x2 = Mod(Q[1],n); y1 = Mod(P[2],n); 
 y2 = Mod(Q[2],n);if(x1==x2&&y1==-y2,return([0]));
 if(P == Q, den=(2* y1 + a[1]* x1 + a[3]); fat=gcd(n,lift(den)); if(fat>1 && fat!=n,return(fat));
 lam = (3*x1^2 + 2*a[2]*x1 + a[4] - a[1]* y1)*den^-1,den= (x2 - x1); fat=gcd(n,lift(den)); if(fat>1 && fat!=n,return(fat)); 
 lam = (y2 - y1)*den^-1); 
 mu = y1 - lam *x1; x3 = lam^2 + a[1]* lam - a[2] - x1 - x2;
 lift([x3, -(lam + a[1])* x3 - mu - a[3]]);
}


{
ellipotn(a,P,m,n)=local(u,t,v);
u = m; t = [0]; v = P;
while(u != 0, if(u%2 == 1, t = ellisumn(a,t,v,n);if(type(t)=="t_INT",return(t))); v = ellisumn(a,v,v,n);if(type(v)=="t_INT",return(v));
u = floor(u/2));
lift(t);
}


{
trovaexp(p,B)=local(ee);
if(p>B,return(0));
ee=1;
while(p^ee<=B,ee++);
ee-1;
}

{
randprimo(ordine)=nextprime(random(10^ordine))
}


{
ecm(n,B)=local(xx,yy,cc,ex,aa,bb,pun,ris);
cc=1;
while(cc==1,
xx=random(n);yy=random(n);aa=random(n);
bb=(yy^2-xx^3-aa*xx)%n;
gg=gcd(4*aa^3+27*bb^2,n);
if(gg>1 && gg!=n, return(gg));
if(gg==1,pun=[xx,yy];cc=0));
forprime(pr=2,B,ex=trovaexp(pr,B);
for(j=1,ex,ris=ellipotn([0,0,0,aa,bb],pun,pr,n);if(ris==[0],return([0]));if(type(ris)=="t_INT",return(ris));pun=ris));
error("fallimento!");
}



{
ellpunti1(a,b,p)=local(e,i,j,punti);
e=ellinit([0,0,0,a,b]);
punti=[];
for(i=0, p-1, for(j=0, p-1, q=[Mod(i,p), Mod(j,p)];
                            if(ellisoncurve(e,q)==1, punti=concat(punti, [q]))));
return(lift(punti))}



{
ellpunti2(a,b,p)=local(e,i,x,s,punti,g,q1,q2);
e=ellinit([0,0,0,a,b]);
punti=[];
for(i=0, p-1, x=Mod(i,p); 
              s=x^3+a*x+b;
              if(kronecker(lift(s),p)==1,g=s^(1/2); q1=[x,g];
                                          q2=[x,-g];
                                          punti=concat(punti,[q1]);
                                          punti=concat(punti,[q2]));
              if(lift(s)==0, q=[x, s]; punti=concat(punti,[q])));
return(lift(punti))}

{
	ellord(z,p)=p+1-ellap(ellinit(z),p)
}	

{fconv(a)=local(n,p,q,pa,pb,qa,qb,lp,lq);n=length(a);pb=1;pa=a[1];lp=[pa];qb=0;qa=1;lq=[qa];for(i=2,n,p=a[i]*pa+pb;q=a[i]*qa+qb;lp=concat(lp,p);
lq=concat(lq,q);pb=pa;pa=p;qb=qa;qa=q);return([lp,lq])}

\\______________________________________________
\\______________________________________________
\\______________________________________________
\\______________________________________________
\\ Programmi fatti per esercizio LEZ 12


\\______________________________________________
\\______________________________________________es 1

\\ Lez 12 Es 1: prodotto di due 
\\ Input: 
\\ Output:
{
prodd(d,a,b)=local(ris);
ris=[a[1]*b[1]+(d^2)*a[2]*b[2],a[1]*b[2]+a[2]*b[1]];
}

\\______________________________________________
\\______________________________________________
\\es 2 funzioni iniziali:

\\ Lez 12 Es 2: da vettore che rappresenta fraz. continua a risultato
\\ Input: V vettore
\\ Output: R frazione corrispondente
{
vett2cf(V)=local(l,R);
l=length(V);
R=V[l];
for(j=1,l-1,R=V[l-j]+1/R);
R;
}

\\ Lez 12 Es 2: da frazione (o eventualmente irrazionale) a frazione continua corrispondente 
\\ Input: F frazione
\\ Output: V vettore corrispondente alla frazione continua di F
{
cf2vett(F)=local();
contfrac(F);
}

\\ Lez 12 Es 2: il singolo p_n della frazione continua V
\\ Input: V
\\ Output: p_n
\\confrontare con conrfracpnqn
{
cfracspn(V,n)=local();
if(n==0, return(V[1]) );
if(n==1, return(V[1]*V[2]+1) );
if(n > 1, return( V[n+1]*cfracspn(V,n-1) + cfracspn(V,n-2)) );
}

\\ Lez 12 Es 2: il vettore dei p_n della frazione continua V 
\\ Input: 
\\ Output:
{
cfracpn(V)=local(l,R);
l=length(V);
R=[];
for(k=0,l-1,R=concat(R,cfracspn(V,k)));
R;
}

\\ Lez 12 Es 2: i singoli q_n della frazioni continue V
\\ Input: 
\\ Output:
{
cfracsqn(V,n)=local();
if(n==0, return(1) );
if(n==1, return(V[2]) );
if(n > 1, return( V[n+1]*cfracsqn(V,n-1) + cfracsqn(V,n-2)) );
}

\\ Lez 12 Es 2: il vettore dei q_n della frazione continua V 
\\ Input: 
\\ Output:
{
cfracqn(V)=local(l,R);
l=length(V);
R=[];
for(k=0,l-1,R=concat(R,cfracsqn(V,k)));
R;
}

\\ Lez 12 Es 2: approxV
\\ Input: frazione continua V come vettore finito.
\\ Output: restituisce la successione pk/qk
{
approxv(V)=local(Vpn,Vqn,R);
Vpn=cfracpn(V); Vqn=cfracqn(V);
R = vector(length(Vpn),k,Vpn[k]/Vqn[k]);
R;
}

\\ Lez 12 Es 2: sviluppo Qn
\\ Input: 
\\ Output:
{
ricQn(d,n)=local();
if(n==0, return(1) );
if(n > 0, return((d-ricPn(d,n-1)^2)/(ricQn(d,n-1))  ) );
}

\\ Lez 12 Es 2: sviluppo Pn
\\ Input: 
\\ Output:
{
ricPn(d,n)=local();
if(n==0, return(0) );
if(n > 0, return(floor((ricPn(d,n-1)+sqrt(d))/(ricQn(d,n-1)))*ricQn(d,n-1)+ricPn(d,n-1)  ) );
}

\\ Lez 12 Es 2: sviluppo x con n
\\ Input: 
\\ Output:
{
xconn(d,n)=local();
return( (ricPn(d,n-1) + sqrt(d))/(ricQn(d,n-1))   );
}

\\ Lez 12 Es 2: sviluppo a con n
\\ Input: 
\\ Output:
{
aconn(d,n)=local();
return( floor(xconn(d,n)));
}

\\ Lez 12 Es 2: quozienti parziali di sqrt d fino a mass
\\ Input: 
\\ Output:
{
seqaconn(d,mass)=local(R);
R=[];
for(i=1,mass,R=concat(R,aconn(d,i)));
R;
}

\\ Lez 12 Es 2: sviluppo Pn
\\ Input: 
\\ Output:
{
ricPn1(d,n)=local();
if(n==0, return(0) );
if(n > 0, return(floor((ricPn(d,n-1)+floor(sqrt(d)))/(ricQn(d,n-1)))*ricQn(d,n-1)+ricPn(d,n-1)  ) );
}

\\ Lez 12 Es 2: sviluppo x con n
\\ Input: 
\\ Output:
{
xconn1(d,n)=local();
return( (ricPn(d,n-1) + floor(sqrt(d)))/(ricQn(d,n-1))   );
}

\\______________________________________________
\\______________________________________________ es 2

\\ Lez 12 Es 2: rad 
\\ Input: riceve d intero non quadrato
\\ Output: sequenza dei quozienti parziali ak di sqrt d con l'algoritmo di Lagrange
\\         lo sviluppo di un irrazionale quadratico e' periodico.  

{
rad(d)=local(a,P,Q,s,v);
s=sqrt(d); v=floor(s);
a=[];P=[0];Q=[1];
k=1;
while(1,a=concat(a,floor((P[k]+v)/Q[k]) );
  P=concat(P,a[k]*Q[k]-P[k]);
  Q=concat(Q,(d-P[k+1]^2)/Q[k]);
  if(a[k]==2*v,break);
  k=k+1);
return(a);

}

\\_______________________________________________
\\ Subordinate a rad:

\\ Lez 12 Es 2: Rad trova la lunghezza del periodo di contfrac(sqrt(d)) 
\\ Input: riceve d intero non quadrato
\\ Output: lungh periodo di contfrac(sqrt(d))  
{
LunghPeriRad(d)=local(V,i,Per);
V=contfrac(sqrt(d));
Per=1;
while(V[Per+1]!=2*V[1],Per=Per+1);
Per;
}

\\ Lez 12 Es 2: Rad trova rad senza usare l'algoritmo di Lagrange, ma usando contfrac 
\\ Input: riceve d intero non quadrato
\\ Output: contfrac(sqrt(d)) con solo 1 periodo 
{
rad2(d)=local(V,i,R);
V=contfrac(sqrt(d));
i=1;
R=[V[1]];
while(V[i]!=2*V[1],R=concat(R,V[i+1]);i=i+1);
R;
}

\\ Lez 12 Es 2: Rad trova rad senza usare l'algoritmo di Lagrange, 
\\              ma usando contfrac e LunghPeriRad 
\\ Input: riceve d intero non quadrato
\\ Output: contfrac(sqrt(d)) con solo 1 periodo 
{
Rad(d)=local();
if(issquare(d)==1,return(contfrac(sqrt(d))) ,
                  return(vector(LunghPeriRad(d)+1,j,contfrac(sqrt(d))[j]) );
);
}

\\ Lez 12 Es 2: Restituisce solo il periodo senza l'antiperiod 
\\ Input: riceve d intero non quadrato
\\ Output: contfrac(sqrt(d)) con solo 1 periodo 
{
PeriRad(d)=local();
vector(LunghPeriRad(d),j,contfrac(sqrt(d))[j+1]);
}

\\ Per avere una verifica è sufficiente confrontare rad con Rad

\\______________________________________________
\\______________________________________________ es 3

\\ Lez 12 Es 3: pell's equation
\\ Input: d
\\ Output:restituisce i valori della soluzione dell'equazione di Pell
\\ devo raddoppiare il periodo: mi serve per le soluzioni dispari
{
pell(d)=local(lungh,X,Y);
lungh=LunghPeriRad(d);
X=cfracpn(concat(rad(d),PeriRad(d)));
Y=cfracqn(concat(rad(d),PeriRad(d)));
if(lungh%2==0,return([X[lungh],Y[lungh]]),
              return([ [X[lungh],Y[lungh]],[X[2*lungh],Y[2*lungh]] ]) )
}

\\________________________________________________
\\________________________________________________ es 4

\\ Lez 12 Es 4: interplus applica n volte la matrice 
\\              [floor(sqrt(d)),d;1,floor(sqrt(d))] al vettore [1;0] 
\\ Input: d intero n numero iteazioni
\\ Output: restituisce la storia dei prodotti.
{
interplus(d,n)=local(R);
R=[;];
\\[floor(sqrt(d)),d;1,floor(sqrt(d))]^n*[1;0];
for(j=1,n, R=concat(R,[floor(sqrt(d)),d;1,floor(sqrt(d))]^j*[1;0]) );
R;
}

\\ Lez 12 Es 4: interplus1 come sopra, ma senza storia, solo il risultato 
\\ Input: 
\\ Output:
{
interplus1(d,n)=local(R);
[floor(sqrt(d)),d;1,floor(sqrt(d))]^n*[1;0];
}

\\ Lez 12 Es 4: interplus come frazioni 
\\ Input: d,n come in interplus
\\ Output: anziché la matrice rest. frazioni.
{
finterplus(d,n)=local(V,R);
V=interplus(d,n); R=[];
for(i=1,n,R=concat(R,V[1,i]/V[2,i]));
R;
}

\\________________________________________________
\\________________________________________________ es 5

\\ Lez 12 Es 5: vettore dei quadrati non perfetti
\\ Input: lim
\\ Output: vettore degli interi non quadrati fino a lim.
{
nonsqu(lim)=local(V,R);
R=[];
V=vector(lim,i,i);
for(j=1,lim,
  if(issquare(j)==0, R=concat(R,j) )
);
R;
}

\\ Lez 12 Es 5: ausiliario condizione (b) Theorem Plus 
\\ Input: h
\\ Output: calcola 
{
ausplus(h)=local();
2*floor(sqrt(h))/(h-floor(sqrt(h))^2);
}

\\ Lez 12 Es 5: determinare i d tali che lo sviluppo ordinario in 
\\              frazioni continue di sqrt(d) abbia periodo minore di 2:
\\              uso la condizione b del teorema Plus
\\ Input: n
\\ Output: lista dei d che soddisfano la condizione (b) fino ad n
{
teoplus(n)=local(Vett,R);
R=[]; Vett=nonsqu(n);
for(i=1,length(Vett),if(type(ausplus(Vett[i]))=="t_INT",R=concat(R,Vett[i])   ));
R;
}

\\________________________________________________
\\________________________________________________ es 6

\\ congettura 1

\\ Lez 12 Es 6: skip 2 pick 3 
\\ Input: V vettore
\\ Output: R vettore con skip 2 pick 3 di V dal primo elemento incluso
\\         torna indietro dal pivot V[5k+1] K=1
{
s2p3i(V)=local(R);
R=[];
R=concat(R,V[1]);
forstep(j=6,length(V),5,R=concat(R,[V[j-2],V[j-1],V[j]]));
R;
}

\\ Lez 12 Es 6: skip 2 pick 3 
\\ Input: V vettore
\\ Output: R vettore con skip 2 pick 3 di V dal primo elemento incluso
\\         va in avanti dal pivot V[5k+1] k=0
{
s2p3a(V)=local(R);
R=[];
R=concat(R,V[1]);
forstep(j=1,length(V)-5,5,R=concat(R,[V[j+3],V[j+4],V[j+5]]));
R;
}

\\ Lez 12 Es 6: skip 2 pick 1, skip 2 pick 3 
\\ Input: V vettore
\\ Output: R vettore con skip 2 pick 1, skip 2 pick 3 di V dal primo elemento incluso
{
s2p1s2p3(V)=local(R);
R=[];
R=concat(R,V[1]);
forstep(j=9,length(V),8,R=concat(R,[V[j-5],V[j-2],V[j-1],V[j]]));
R;
}

\\ Si in avanti che indietro funzionano correttamente, ma 
\\ solo con vettori sufficientemete lunghi. Aumento i periodi
\\questo però appesantisce troppo i calcoli 
\\in caso di periodi troppo lunghi. Sarebbe da migliorare.

\\ Lez 12 Es 6: approssimazione pn/qn per rad(d)
\\ Input: d
\\ Output: pn/qn per rad(d) con 3 periodi
{
approfc(d)=local(P,Q,dlungo,R);
dlungo=concat(rad(d),concat(PeriRad(d), PeriRad(d)));
P=cfracpn(dlungo);
Q=cfracqn(dlungo);
R=vector(length(P));
for(j=1,length(P),R[j]=P[j]/Q[j]);
R;
}

\\ Lez 12 Es 6: confronto la frazione continua e la matrice convergente  
\\ Input: d intero squarefree
\\ Output: convergenza frazione continua e matrice convergente
\\         uso approxv e finterplus di sqrt(d)
{
confronto(d)=local(A,B);
A=approfc(d);
B=finterplus(d,length(A));
return([A;B]);
}

\\ Lez 12 Es 6: confronto applicando s2p3i o s2p1s2p3  
\\ Input: d intero squarefree che soddisfi alla condizione di conjecture 1
\\ Output: valuta i due vettori della congettura
{
congett2(d)=local(C,S);
C=confronto(d);
if(C[1,1][1]%2==0,S=s2p1s2p3(C[1,1]),S=s2p3i(C[1,1]));
return([C;S]);
}

\\________________________________________________
\\________________________________________________ es 7

\\ Lez 12 Es 7: modifica di ceil
\\ Input: a reale
\\ Output: floor(a)+1 come da definizione pag 1
{
ceilm(a)=local();
floor(a)+1;
}


\\ Lez 12 Es 7: tassello dello sviluppo f.c. minus 
\\ Input: a
\\ Output: 1/(ceilm(a)-a)
{
gaussm(a)=local();
1/(ceilm(a)-a);
}

\\ Lez 12 Es 7: i primi n termini dello sviluppo di a in f c minus 
\\ Input: a numero da sviluppare in f c, n numero di termini in f c
\\ Output: fc minus di a con n termini
\\         uso Beta come vettore parallelo seguendo pag 4 senza ricorsione.
{
contfracm(a,n)=local(Beta,B);
Beta=vector(n);
B=vector(n);
B[1]=ceilm(a);

Beta[1]=a;
Beta[2]=1/(B[1]-Beta[1]);

for(i=2,n,Beta[i]=gaussm(Beta[i-1]); B[i]=ceilm(Beta[i]) );
B;
}



\\________________________________________________
\\________________________________________________ es 8


\\ Lez 12 Es 8: calcolo degli Rn di minus 
\\ Input: V vettore
\\ Output: vettore degli Rn
{
fcR(V)=local(l,R);
l=length(V)+1; R=vector(l);
R[1]=0;
R[2]=1;
for(i=3,l,R[i]=V[i-2]*R[i-1]+R[i-2]);
R;
}

\\ Lez 12 Es 8: calcolo degli Sn di minus 
\\ Input: V vettore
\\ Output: vettore degli Sn
{
fcS(V)=local(l,S);
l=length(V)+1; S=vector(l);
S[1]=0;
S[2]=1;
for(i=3,l,S[i]=V[i-2]*S[i-1]+S[i-2]);
S;
}

\\ Lez 12 Es 8: 
\\ Input: V sviluppo minus di un valore a
\\ Output: relativi Rk ed Sk dello sviluppo
{
fconvm(V)=local();
[fcR(V);fcS(V)];
}


\\RICORSIONE:
{
ricmR(v,n)=local();
if(n == 0,return(0));
if(n == 1, return(1));
if(n > 1, return(V[n-1]*ricmR(V,n-1)-ricmR(V,n-2)  ));
}

{
ricmS(v,n)=local();
if(n == 0,return(-1));
if(n == 1, return(0));
if(n > 1, return(V[n-1]*ricmS(V,n-1)-ricmS(V,n-2)  ));
}



\\________________________________________________
\\________________________________________________ es 9

\\ Lez 12 Es 5: ausiliario condizione (b) Theorem Minus 
\\ Input: h
\\ Output: calcola 
{
ausminus(h)=local();
2*ceilm(sqrt(h))/(ceilm(sqrt(h))^2-h);
}

\\ Lez 12 Es 5: determinare i d tali che lo sviluppo ordinario in 
\\              frazioni continue di sqrt(d) abbia periodo minore di 2:
\\              uso la condizione b del teorema Plus
\\ Input: n intero
\\ Output: lista dei d che soddisfano la condizione (b) fino ad n
{
teominus(n)=local(Vett,R);
R=[]; Vett=nonsqu(n);
for(i=1,length(Vett),if(type(ausminus(Vett[i]))=="t_INT",R=concat(R,Vett[i])   ));
R;
}



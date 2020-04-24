
/*

LABORATORIO DI CALCOLO SIMBOLICO
ANNO ACCADEMICO 2010/2011

SEBASTIANO FERRARIS

\r Ferraris\lezione_9

*/

\\____NONA LEZIONE_________________________________________________
\\____________________________________________________________________
\\____________________________________________________________________

\\____________________________________________________________________
\\Programmi dati in classe LEZ 9

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
for(j=1,n,resti[j]=[];
  for(i=1,s,resti[j]=concat(resti[j],pol2vetc(divrem(x^(j-1),poli[i])[2],poldegree(poli[i])))));
for(i=1,n,T[,i]=mattranspose(resti[i]));
T;
}

{
\\applica del(n) come pura trasformazione lineare ad un vettore
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


\\_____________________________________________
\\ Programmi fatti per esercizio LEZ 9

\\________________________________________________________________
\\________________________________________________________________es 1

\\ Lez 9 Es 1: ausiliaria di gam, prende un vettore di tau vettori, 
\\            un vettore di interi positivi lungo tau e forma la matrice le cui
\\            colonne sono i vettori del vettore dato ripetuto nella matrice e
\\            shiftato tante volte quante lo dice il vettore dato  
\\ Input: vett, grad
\\ Output: matrice M t.c. M=destrak(vettore[j],grad[j+h]) 
{
incoleshift(vett,n,grad,tau)=local(M,VN);
M=matrix(n,n);
VN=[];
for(h=1,tau, for(k=1,grad[h],VN=concat(VN,destrak(vett[h],k-1))  ));
for(i=0,n-1, for(j=1,n, M[j,i+1]=VN[j+i*n]));
M
}


\\ Lez 9 Es 1: trova la rappresentazione della funzione gammaQ come matrice
\\ Input: n esponente di x^n-1
\\ Output: restituisce la matrice G che rappresenta la trasformazione lineare
\\         data da gammaQ
{
gam(n)=local(G,fatt,tau,vgra);

fatt=factor(x^n-1);
tau=matsize(fatt)[1];
G=vector(tau);

\\vgra è il vettore dei gradi lungo tau:
vgra=[];
for(i=1,tau,vgra=concat(vgra,poldegree(fatt[i,1]) ));

\\vettore G: è un vettore lungo tau i cui elementi 
\\sono i vettori colonna di G senza shift  

for(i=1,tau,G[i]= pol2vetc(gammaQ(matid(tau)[i,],n),n) );

\\applico incoleshift a G ed ottengo la matrice cercata!
incoleshift(G,n,vgra,tau);
}

\\ Lez 9 Es 1: applica gam(n) come pura trasformazione lineare ad un vettore
\\ Input: vettore a del prodotto dei campi quozientati sui poly ciclotomici
\\ Output: applica la funzione gamma ad a tramite la matrice costruita sopra
{
gamlin(a)=local(n);
n=length(a);
mattranspose(gam(n)*mattranspose(a));
}

\\_____________________________________________________________
\\_________________________________________________es 2 parte 1
\\ Lez 9 Es 2: ausliaria 
\\ Input: vettore di polinomi, lunghezza ed un vettore di gradi 
\\ Output: trasforma la sequenza dei polinomi in un vettore dei coefficienti
\\         seguendo i gradi dati dal vettore dei gradi 
{
polinomizza(vett,lungh,vettgrad)=local(ris);
ris=[];
for(i=1,lungh,ris=concat(ris,pol2vetc(vett[i],vettgrad[i])));
ris;
}

\\ Lez 9 Es 2: delnp(n,p) restituisce la trasformazione lineare di delta
\\ Input: n p 
\\ Output: matrice che rappresenta la trasformazione lineare
\\         data da delnp
{
delnp(n,p)=local(fatt,tau,vgra,T);
if(n==p,T=Mod(matid(n),n),

fatt=factormod(x^n-1,p);
tau=matsize(fatt)[1];

\\vettore dei gradi:
vgra=[];
for(i=1,tau,vgra=concat(vgra,poldegree(fatt[i,1]) ));

\\matrice T
T=matrix(n,n,i,j, Mod(  polinomizza(teta(Mod(1,p)*x^(j-1),n,p),tau,vgra)[i]  ,p)  );
);

T;
}

\\ Lez 9 Es 1: applica delnp(n,p) come pura trasformazione lineare ad un vettore
\\ Input: vettore a di un elemento Z_p[x]/(x^n-1)
\\ Output: applica la funzione delnp ad a tramite la matrice costruita sopra
{
delnplin(a,p)=local(n);
n=length(a);
mattranspose(delnp(n,p)*mattranspose(a));
}

\\Nota: polinomizza per n=p non funziona, motivo per cui ho messo if in delnp che
\\eviti quel caso. Sarebbe più elegante usare un algoritmo che non abbia la funzione
\\ausiliaria polinomizza.

\\_____________________________________________________________
\\_________________________________________________es 2 parte 2

\\ Lez 9 Es 2: trova la rappresentazione della funzione gamm come matrice
\\ Input: n esponente di x^n-1, p modulo
\\ Output: restituisce la matrice G che rappresenta la trasformazione lineare
\\         data da gamm
{
gamnp(n,p)=local(R,G,fatt,tau,vgra);
if(n==p,R=Mod(matid(n),n),

fatt=factormod(x^n-1,p);
tau=matsize(fatt)[1];
G=vector(tau);

vgra=[];
for(i=1,tau,vgra=concat(vgra,poldegree(fatt[i,1]) ));
 
for(i=1,tau,G[i]= Mod(pol2vetc(gamm(matid(tau)[i,],n,p),n) ,p)   );

R=incoleshift(G,n,vgra,tau);
);
R
}

\\ Lez 9 Es 1: applica gamnp(n) come pura trasformazione lineare ad un vettore
\\ Input: vettore a del prodotto dei campi quozientati sui poly ciclotomici
\\ Output: applica la funzione gamma ad a tramite la matrice costruita sopra
{
gamnplin(a,p)=local(n);
n=length(a);
mattranspose(gamnp(n,p)*mattranspose(a));
}

\\Nota: incoleshift per n=p non funziona, motivo per cui ho messo if in delnp che
\\eviti quel caso. Sarebbe più elegante usare un algoritmo che non abbia la funzione
\\ausiliaria incoleshift.

\\____________________________________________________________
\\______________________________________________________Es 3


\\ Lez 9 Es 3: ausiliaria, dato un numero dec ed una lunghezza lungh 
\\             restituisce dec in binario in un vettore binario lungo lungh  
\\ Input: dec, lungh
\\ Output: rv vett, binario
{
binarylong(dec,lungh)=local(bv,rv);
rv=vector(lungh);
bv=binary(dec);
l=length(bv);
for(i=1,l,rv[i]=bv[l-i+1]);
rv;
}

\\ Lez 9 Es 3: funzione hammH che retituisce la matrice di controllo del codice
\\             lineare di Hamming, di altezza m. (Non è in forma standard)
\\ Input: m
\\ Output: Matrice h di controllo
{
hammH1(r)=local(n,v,M);
n=2^r-1;
v=vector(n);
for(i=1,n,v[i]=binarylong(i,r));
M=matrix(r,n,i,j,v[j][i]);
M;
}

\\Cerco un sistema per trovare la matrice H in forma standard
\\costruendo un vettore del tipo [1,2,4,...,2^k,3,5,6,7...,n]:

\\ Lez 9 Es 3: ausilairia, elimina l'elemento x dal vettore v 
\\ Input: v vettore, x elemento (necessariamente) in v
\\ Output: v senza x
{
elimvx(v,x)=local(h,pos,vris);
h=length(v);
pos=1;
i=1;
vris=vector(h-1);
while(v[i]!=x,i=i+1;pos=i);
if(pos==1, for(l=1,h-1,vris[l]=v[l+1]),
   for(j=1,pos-1,vris[j]=v[j]);
   for(k=pos+1,h,vris[k-1]=v[k]) );
vris;
}

\\ Lez 9 Es 3: ausliaria, elimina dal vettore v gli elementi del vettore E 
\\ Input: v vettore, eE vettore di elementi contenuti in v
\\ Output:
{
elimvxvett(v,E)=local(h,ris);
h=length(E);
ris=v;
for(j=1,h,ris=elimvx(ris,E[j]));
ris;
}

\\ Lez 9 Es 3: genera il vettore [1,2,4,...,2^k,3,5,6,7...,n]
\\ Input: n
\\ Output: [1,2,4,...,2^k,3,5,6,7...,n] k=floor(log2(n))
{
vetthamm(n)=local(v1,v2);
v1=vector(n,i,i);
v2=[];
for(j=0,floor(log2(n)),v2=concat(v2,2^j));
v1=elimvxvett(v1,v2);
concat(v1,v2);
}


\\ Lez 9 Es 3: funzione hammH che retituisce la matrice di controllo del codice
\\             lineare di Hamming, di altezza m, in forma standard!
\\ Input: m
\\ Output: Matrice h di controllo
{
hammH(r)=local(n,v,vh,M);
n=2^r-1;
v=vector(n);
vh=vetthamm(n);
for(i=1,n,v[i]=binarylong(vh[i],r));
M=matrix(r,n,i,j,v[j][i]);
M;
}

\\Esiste un sistema per ricavare la stessa matrice dalla funzione teta, o in un 
\\modo meno arzigogolato?

\\_________________________________________________________________
\\_________________________________________________________________Es 4

\\ Lez 9 Es 4: hammG che restituisce la matrice G in forma standard del codice
\\             lineare di Hamming
\\ Input: r altezza del codice
\\ Output: matrice generatrice del (2^r-1,2^r-1-r)-codice di Hamming
{
hammG(r)=local(mid,H,ris);
ris=matrix(2^r-1-r,2^r-1);
mid=matid(2^r-1-r);
H=hammH(r);
for(i=1,2^r-1-r,ris[i,]=concat(mid[i,],H[,i]~));
ris;
}

\\______________________________________________________________
\\______________________________________________________________Es5

\\ Lez 9 Es 5: codifica di un vettore con il codice di Hamming	
\\ Input: r altezza del codice, a vettore (di lunghezza k=2^r-1-r)
\\ Output: vettore codificato
{
hcode(r,a)=local();
lift(Mod(a*hammG(r),2))
}

\\_______________________________________________________________
\\_______________________________________________________________Es 6

\\lez 9 es 6: ausiliaria, dato un vettore ridà solo i primi k elementi dello stesso,
\\Input: v vettore, k <= length v
\\Output: ris vettore avente i primi k elementi di v
{
prendik(v,k)=local(ris);
ris=[];
for(i=1,k,ris=concat(ris,v[i]));
ris;
}

\\lez 9 es 6: ausiliaria, cambia un bit del vettore alla posizione ps
\\Input: V vettore binario s indice di un elemento del vettore
\\Output: nV vettore binario con elemento posto j invertito
{
correggipos(V,s)=local();
V[s]=lift(Mod(V[s]+1,2));
V;
}

\\Lez 9 es 6: ausilairio, sindrome di un vettore 
\\(usata solo come prova e non in hdecoder, altrimenti calcola due volte H)
{
hsind(r,x)=local(H,sind);
H=hammH(r);
sind=lift(Mod(H*x~,2));
sind;
}

\\ Lez 9 Es 6: decodifica un vettore con il codice di hamming
\\ Input: r altezza del codice, a vettore in codice (di lunghezza n=2^r-1) 
\\ Output: vettore decodificato
{
hdecode(r,x)=local(H,sind,k,pos);
k=2^r-1-r;
H=hammH(r);
sind=lift(Mod(H*x~,2));
pos=0;
if(sind==vector(r)~, return(prendik(x,k)), 
    for(j=1,matsize(H)[2],if(sind==H[,j],pos=j)              
    )
);
if(pos==0,
    return("due o + errori"),
   return(prendik(correggipos(x,pos),k)) );
}

/*
Commento sull'algortimo:
confronta la sindrome con i vettori colonna della matrice h. 
Se la sindrome è il vettore colonna nullo, allora non si presume non ci siano errori.
Se la sindrome non è il vettore nullo, allora la si confronta con le 
prime k colonne della matrice h.
Se la sindrome è la j-esima colonna allora si corregge risp alla posizione j
Se la sindrome è diversa da tutte le colonne, allora il codice non corregge.
Esempio del suo uso

(13:40) gp > v=[1,0,1,0]
%4 = [1, 0, 1, 0]

(13:40) gp > vcod=hcode(3,v)
%5 = [1, 0, 1, 0, 1, 0, 1]

(13:40) gp > verr=[1,0,1,0,0,0,1]
%6 = [1, 0, 1, 0, 0, 0, 1]

(13:40) gp > hdecode(3,vcod)
%7 = [1, 0, 1, 0]

(13:40) gp > hdecode(3,verr)
%8 = [1, 0, 1, 0]

*/






/*

LABORATORIO DI CALCOLO SIMBOLICO
ANNO ACCADEMICO 2010/2011

SEBASTIANO FERRARIS

\r Ferraris\lezione_8

*/

\\____Ottava LEZIONE_________________________________________________
\\____________________________________________________________________
\\____________________________________________________________________

\\____________________________________________________________________
\\Programmi dati in classe LEZ 8


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

{\\schema dell'algoritmo potenza veloce; operid � l'elemento neutro, oper � l'operazione
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

{\\restituisce 1 se n � il periodo di elem (con operazione oper e neutro operid, vedi pot)
periodo(elem,n)=local(ff,s);
if(pot(elem,n)!=operid,return(0));
ff=factor(n);
s=matsize(ff)[1];
for(i=1,s,if(pot(elem,n/ff[i,1])==operid,return(0)));
1;
}

{\\trova il periodo di elem (con operazione oper e neutro operid, vedi pot) se ordine � un suo multiplo
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

{\\riceve due polinomi a,b in x, di grado <= n-1 e calcola a*b mod x^n-1 , dove n � un intero positivo
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

{\\inversa di Fourier; fournp(n,p)*Mod(circol(vettore),p)*fournpinv(n,p) � diagonale
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

{\\d� ric(a,b,h,k,n) modulo m, utilizzando la matrice compagna
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
\\riceve tre polinomi a,b,p in una stessa variabile e calcola a*b mod (p,n) , dove n � un intero positivo
moltp(a,b,p,n)=local(y,z);
y=Mod(Mod(1,n)*a,p);z=Mod(Mod(1,n)*b,p);
return(lift(lift(y*z)))}

{
\\riceve due polinomi in una stessa variabile e due interi e=esponente, n=modulo
\\n deve essere positivo, e pu� essere negativo, se a � coprimo con p
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
\\restituisce un array ci s polinomi, dove s � il numero dei poli irriducibili nei quali x^n-1 si fattorizza modulo p
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
decmatprodQ(ma)=local(ee,s,a,n,vv);
a=Polrev(ma[1,]);n=matsize(ma)[1];
ee=circolidemQ(n);s=length(ee);
vv=vector(s,k,1+moltnQ(ee[k],a-1,n)); \\he per def la comp moltiplicativa
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
\\probabilit� che una circolante nxn abbia il determinante divisibile per il primo p, n generico
npprobgen(n,p)=local();
if(n%p==0,return(0));
1.-circolgrorder(n,p)/(p^n);
}

{
\\probabilit� che una circolante nxn abbia il determinante divisibile per il primo p, posto che n divida p-1
npprob(n,p)=local();
if((p-1)%n != 0,return(0));
1.-(1-1/p)^n;
}

{
\\probabilit� ricavata sperimentalmente
npprobex(n,p,casi)=local();
favorevoli=0;
for(i=1,casi,if(matdet(circol(randvett(n,p)))%p==0,favorevoli++));
1.*favorevoli/casi;
}


\\_____________________________________________
\\_____________________________________________
\\ Programmi fatti per esercizio LEZ 8

\\_____________________________________________
\\ Lez 8 Es 1: Test di Miller
\\ Input: n, a
\\ Output: test di miller di n sulla base a
{
miller(n,a)=local(mcd,risp);
mcd=gcd(n,a);
if(mcd!=1, ris=0, if( lift(Mod(a^(n-1),n))==1, ris=1,ris=0));
ris;
}

\\_____________________________________________
\\ Lez 8 Es 2: Test di Baillie
\\ Input: n intero
\\ Output: esegue il test di baillie su n. Consiste nel test di Miller in base 2
\\         seguito dal test di Lucas
{
baillie(n)=local(mill,ris);
mill=miller(n,2);
if(mill==0,ris=0,ris=lucasB(n)); \\oppure uso lucas(n,1000)
ris;
}

{
lucasB(n)=local(D,jj,U,risp);
D=5;
jj=1;
while(kronecker(D,n)!=-1,D=((-1)^(jj))*(abs(D)+2);jj=jj+1  );
U=ric(0,1,1,(D-1)/4,n+1);
if(U%n==0,ris=1,ris=0);
ris;
}

\\______________________________________________
\\ Lez 8 Es 3: aks passo 2
\\ Input: n intero
\\ Output: si calcola l'algoritmo aks fino al secondo passo
{
akspasso2(n)=local(r,mcd,logar);
if(issquare(n)==1,return(0),
   r=2;
   while(r<n,if(gcd(n,r)>1,return(0),
              for(i=1,4*floor(log(n)^2),if(Mod(n^r,r)!=1,break ) )
              );
   r=r++;
   );

);
\\return(r); \\debug
}

\\______________________________________________
\\ Lez  Es : aks
\\ Input: n intero
\\ Output: aks per n
{
aks(n)=local(r);
if(issquare(n)==1,return(0),
   r=2;
   while(r<n,if(gcd(n,r)>1,return(0),
              for(i=1,4*floor(log(n)^2),if(Mod(n^r,r)!=1,break ) )
              );
   r=r++;
   );

);

if(r==n, return(1),
   for(j=1,2*floor(sqrt(r))*floor(log(n))+1,
        if(Mod( (Mod(1,n)*x+Mod(j,n))^n, x^r-1)== Mod(Mod(1,n)*x^n+Mod(j,n)^n,x^r-1),
            return(0), return(1)
        );
   );

);
}

\\"lento" ma funziona:
\\ 7 secondi per 2333
\\ >56 sec per 20719


\\______________________________________________
\\ Lez 8 Es 5: del(n)
\\ Input: n intero
\\ Output: restituisce la matrice T che rappresenta la trasformazione lineare
\\         data da TetaQ
{
del(n)=local(fatt,tau,vgra,T);
fatt=factor(x^n-1);
tau=matsize(fatt)[1];

\\vettore dei gradi:
vgra=[];
for(i=1,tau,vgra=concat(vgra,poldegree(fatt[i,1]) ));
\\vgra; \\vgra per debug

\\matrice T
T=matrix(n,n,i,j,polinomizza(tetaQ(x^(j-1),n),tau,vgra)[i]);
T;
}


{
\\rest la iesima colonna della matrice T, con
\\    vett=tetaQ(x^i,n)
\\    lungh=tau
\\    vettgrad=vgra vettore dei gradi dei poly della fattor di x^n-1

polinomizza(vett,lungh,vettgrad)=local(ris);
ris=[];
for(i=1,lungh,ris=concat(ris,pol2vetc(vett[i],vettgrad[i])));
ris;
}


























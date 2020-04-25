/*

LABORATORIO DI CALCOLO SIMBOLICO
ANNO ACCADEMICO 2010/2011

Gli esercizi ridondandti e le prove sono commentati.

SEBASTIANO FERRARIS

\r Ferraris\lezione_5

*/

\\____Quinta LEZIONE_________________________________________________
\\____________________________________________________________________
\\____________________________________________________________________

\\____________________________________________________________________
\\Programmi dati in classe LEZ 5



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
gaps1(a,b,m)=local(lis);
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
write(nomefile,a)
}

\\ Lez 1 Es2 ter: Euclide Mullin. Stampa dei passaggi fino a niter.
\\ Input: riceve vettore (a_i) di interi e il numero di iterazioni desiderate.
\\ Output: no. Stampa i vettori concatenati con il fattore minimo primo di Y = prod(a_i) + 1
{
euc1seq(V,niter) =  local();
print(V);
for(m = 0, niter - 1,  V = concat(V,factor(prod(n = 1, length(V), V[n]) +1)[1,1]  ); print(V) );
}

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
if(e>0,b=floor(n^(1/e));if(b==2 || !isprime(b), return([]));
for(i=2,n-2,if(gcd(i,n)==1 && znorder(Mod(i,n))==phi,ris=concat(ris,i)));
return(ris));
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

{\\dati due vettori a,b restituisce il vettore (ordinato e senza ripetizioni) formato dagli elementi di a che non sono in b
complemento(a,b)=local();
vecsort(eval(setminus(Set(a),Set(b))));
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

{\\calcola l'insieme delle orbite di Z_n sotto l'azione del grupo contenuto in (Z_n)*
classi(n,gruppo)=local(nn,X,G,ris,gg,ele,orbita);
X=vector(n,k,k-1);gg=length(gruppo);if(gg==0,return(0),G=gruppo);ris=[];nn=n;
while(nn!=0,ele=Mod(X[1],n);orbita=[];
for(ii=1,gg,orbita=concat(orbita,G[ii]*ele));
orbita=elimina(lift(orbita));ris=concat(ris,[orbita]);X=complemento(X,orbita);nn=length(X));
ris;
}

\\_______________________________________________________________
\\_______________________________________________________________
\\ Programmi fatti per esercizio LEZ 5 (i primi 3 es della scheda)

\\
\\_______________________________________________________________
\\ lezione 5 esercizio 1
\\

\\ Lez 5 Es 1 parte 1: trova h tale che kronecker(h^2 - 4*a,p)=-1
\\ ...si va per tentativi...
\\ Input: a, p tali che kronecker(a,p) = 1
\\ Output: restituisce l'h tale che kronecker(h^2 - 4*a,p)= -1
{
trovah(a,p) = local(h);
h = 4*a+1;
while(kronecker(h^2 - 4*a,p) != -1, h=h+1);
h;
}

\\ Lez 5 Es 1 parte 2: qsolve(a,p) risolve x^2 = a (mod p) con il metodo di lucas
\\ Input: a e p primo, da x^2 = a (mod p)
\\ Output: se esistono le soluzioni della congruenza
\\ W_n(a,b,h,k)=ric(a,b,h,k,n)
\\ V_n(h,k) = W_n(2,h,h,k) = ric(a,b,h,k,n)

{
qsolve(a,p) = local(hh, sol1, vettsol);
hh = trovah(a,p);
if(kronecker(a,p)!=1 ,print("No soluzioni"), sol1 = (1/2)*ric(2,hh,hh,a,(p+1)/2);
                                             vettsol=[Mod(sol1,p),- Mod(sol1,p)];
                                             vettsol;
);
}

\\
\\_______________________________________________________________
\\ lezione 5 esercizio 2
\\

\\ Lez 5 Es 2: parte 1 fibonacci generalizzata ricorsione
\\ Input: (a,b,h,k,n)
\\ Output: W_n(a,b,h,k)
{
fibgen(a,b,h,k,n) = local();
if(n == 0, return(a));
if(n == 1, return(b));
if(n > 1, return( h*fibgen(a,b,h,k,n-1) - k*fibgen(a,b,h,k,n-2)) );
}

\\ Lez 5 Es 2: parte 2 fibonacci generalizzata vettore
\\ Input: (a,b,h,k,n)
\\ Output: il vettore di tutti gli elementi W_i(a,b,h,k) i <= n
{
fibgenvett(a,b,h,k,n) = local(vett,cont);
vett=[];
for(cont = 0, n, vett = concat(vett,fibgen(a,b,h,k,cont)));
vett;
}

\\ Lez 5 Es 2: parte 2 fibonacci generalizzata con a=0,b=1, cioè U_n (vettore)
\\ Input: (h,k,n)
\\ Output: il vettore di tutti gli elementi U_i(a,b,h,k) i <= n
{
fibgenvett_U(h,k,n)=local(risp);
fibgenvett(0,1,h,k,n);
}

\\ Lez 5 Es 2: parte 3 apparizione che restituisce il rango di apparizione di n in U_m(h,k)
\\ Input: h,k,n che permette di calcolare la sequenza U_i(h,k), i<=m
\\ Output: rango di apparizione di n, cioè il più piccolo m tale che
\\         n divide U_i(h,k)
{
apparizione(h,k,n)=local(Uvett,m);
m=1;
while(fibgen(0,1,h,k,m)%n != 0, m = m+1);
m;
}

\\NB troppo lento!

\\
\\_______________________________________________________________
\\ lezione 5 esercizio 3
/*
calcolare alcuni elementi della successione apparizione(1,-1,p^n), con n = 1,2,...,s
e con p primo dispari, e formulare una congettura!
*/
\\ Lez 5 Es 3: congettura sul rango di apparizione di p^n in U_m(h,k)
\\ Input: (h,k,p,s)
\\ Output: calcola apparizione(1,-1,p^n), con n = 1,2,...,s e con p primo dispari

\\ Tempo di calcolo di apparizione(1,-1,3^3) maggiore di 30 sec.
{
congettapparizione(h,k,p,s)=local(vett);
vett=[];
for(i=1,s,vett = concat(vett,apparizione(h,k,p^i)));
vett;
}

/*
Per poter formulare una congettura serve un algorimo più efficiente...
TENTATIVO:

[  n|F_m <=> n|m && m|F_m ]
condizone da scoprire è quando m|F_m
*/
\\ Input: m
\\ Output: se m|F_m ritorna 1 altrimenti 0
{
divfib(m) = local();
if(fibgen(0,1,1,-1,m)%m == 0, return(1), return(0));
}

\\ Input: n e mettere m=1 (ricorsione)!
\\ Output: m tale che n|F_m
{
apparizione1(n,m)=local();
if(divfib(m)==1 && m%n == 0, return(m), apparizione1(n,m+1));
}

{
congettapparizione1(p,s)=local(vett);
vett=[];
for(i=1,s,vett = concat(vett,apparizione1(p^i,1)));
vett;
}

/*
ancora lento dato che divfib è lento come apparizione!
cercare un sistema per trasformare F_m in una somma, e vedere se m divide i singoli
addendi...
*/


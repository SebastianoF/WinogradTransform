/*

LABORATORIO DI CALCOLO SIMBOLICO
ANNO ACCADEMICO 2010/2011

SEBASTIANO FERRARIS

\r Ferraris\esercizi

*/


\\____QUARTA LEZIONE_________________________________________________
\\____________________________________________________________________
\\____________________________________________________________________

\\____________________________________________________________________
\\Programmi dati in classe LEZ 4




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

{
periodo(elem,n)=local(ff,s);
if(pot(elem,n)!=operid,return(0));
ff=factor(n);
s=matsize(ff)[1];
for(i=1,s,if(pot(elem,n/ff[i,1])==operid,return(0)));
1;
}

{
cercaperiodo(elem,ordine)=local();
fordiv(ordine,n,if(pot(elem,n)==1,return(n)));
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


\\Esercizi NUOVI della quarta lezione

{
mersenne(p)=local(n,a);
n=2^p-1;a=Mod(4,n);
for(i=1,p-2,a=a^2-2);
if(lift(a)==0,return("OK"),return("NO"))}

\\test di lucas nella prima forma, cerca anche il limite di a, base
\\massima del test di Lucas
{
lucas(n,lim)=local(qf,zz,le,ris);
qf=factor(n-1)[,1];zz=vector(lim,k,prime(k));
le=length(qf);
for(i=1,le,qq=qf[i];ris=lucascond(n,qq,zz);
if(ris==0,print(qq);return(0));
if(ris==-1,return(-1)));
return(1)
}

\\verifica se è soddisfatta la condizione di lucas per i tre valori dati
\\
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

\\ fabbrica n=i*F+1 insiemi di primi molto specializzati!
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

\\lavora in Z_n[x]/(x^r -1) , (x+a)^n = x^n +a^n
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

\\Fa proprio a*b con il prodotto delle permutazioni

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
permatinv(a)=local(ris,n);
n = matsize(a)[1];
ris=vector(n);
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

\\_____________________________________________
\\ Programmi fatti per esercizio LEZ 4


\\ Lez 4 Es 1: permper(a) riceve una permutazione e restituisce il suo periodo.
\\ usare il formato delle permutazioni descritto in percomp 
\\ Input: permutazione a
\\ Output: periodo di a

\\tentativo con cercaperido
/*
{
permper(a)=local(b,c);
oper(b,c)=permcomp(b,c);
operid=permid(length(a));
cercaperiodo(a,length(a));
}
*/

\\senza cercaperiodo:
{
permper(a)=local(lungh,cont,b);
b=a;
cont = 1;
lungh =length(a);
while(b != permid(lungh), b = permcomp(a,b); cont=cont+1);
cont;
}

\\ Lez 4 Es 2: mat2pol dalla matrice circolante associa il polinomio isomorfo 
\\ Input: matrice circolante
\\ Output: polinomio ismomrfo in F[x]/(x^n - 1)
{
mat2pol(a)=local(vett,ris);
vett = a[1,];
ris = Mod(Polrev(vett),x^(length(vett))-1);
ris; 
}

\\ Lez 4 Es 2 bis: pol2mat(b,n)
\\ Input: polinomio b da mettere in F[x]/(x^n - 1), con n dall'utente
\\ Output: matrice circolante corrispondente
{
pol2mat(b,n)=local(modpol,vett,ris);
modpol = Mod(b, x^n - 1);
vett=pol2vect(lift(modpol),n);  
ris=circol(vett);
ris;
}


\\ Lez 4 Es 3: matinvQ(a)
\\ Input: matrice circolante ad elementi razionali
\\ Output: 0 se non è inv 1 se è inv, senza calcolare il det (si usa l'isomorfismo)
\\ la somma di ogni termine del polinomio è uguale al prodotto dei termini per il grado,
\\ allora il det è uguale a zero.
 
{
matinvQ(a)=local(vett,n,hh,kk);
vett=a[1,];
n=length(vett);
hh=sum(i=1,n,vett[i]^2);
kk=n*prod(i=1,n,vett[i]);
if(hh==kk,return(0),return(1));
}

\\ Da migliorare capendo dove il polinomio al quadrato calcolato in 1 è uguale
\\ alla somma dei termini al quadrato. 


\\ Lez 4 Es 4: radnp(n,p) restituisce una radice primitiva n-esima di 1 modulo p
\\ Input: n, p
\\ Output: resituisce quel numero alpha che elevato ad n è radice primitiva 
\\ modulo p. Cioè ord_p(alpha^n) = phi(p) = p-1
{
radnp(n,p) = local(zz,aa);
zz = znprimroot(p);
if((p-1)%n == 0, aa = zz^((p-1)/n); return(aa), return("errore"));
}

\\ Lez 4 Es 5:  fournp(n,p). Sia alpha_n = radnp(n,p). restituisce la matrice
\\avente come elemento al posto i,k (alpha_n)^(i*j)
\\ Input: n intero, p primo opportuni (?)
\\ Output: matrix(n,n,i,j,(alpha_n)^((i-1)*(j-1)))
{
fournp(n,p) = local(alpha, MM);
alpha = radnp(n,p);

\\solo il vettore di tutte le soluzioni:
\\MM=vector(n,k,alpha^k);

\\la matrice come richiesta dall'esercizio
MM = matrix(n,n,i,j,(alpha)^((i-1)*(j-1)));
MM;
}

\\ Lez 4 Es 5 Bis: fournpinv(n,p) che calcola la matrice inversa di fournp(np)
\\ Input: n,p
\\ Output: calcola MM~ dove MM = fournp(n,p)
{
fournpinv(n,p) = local(MM, MMtrasp);
MM = fournp(n,p);
MMtrasp = MM~;
MMtrasp;
}


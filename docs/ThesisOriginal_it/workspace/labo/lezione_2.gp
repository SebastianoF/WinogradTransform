/*

LABORATORIO DI CALCOLO SIMBOLICO
ANNO ACCADEMICO 2010/2011

SEBASTIANO FERRARIS

\r Ferraris\esercizi

*/



\\____SECONDA LEZIONE_________________________________________________
\\____________________________________________________________________
\\____________________________________________________________________

\\____________________________________________________________________
\\Programmi dati in classe LEZ 2

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

\\importante per salvare il risutlato in un file .txt !! comando write nomefile: zzz.gp
{\\usa euc1 n volte e scrive il risultato nel file nomefile ES: useuc1([2],16,"cerruti/eu16")
useuc1(lista,n,nomefile)=local(a);
a=lista;
for(i=1,n,a=euc1(a));
write(nomefile,a)
}

{\\schema dell'algoritmo potenza veloce; operid è l'elemento neutro, oper è l'operazione
\\NB oper e operid vanno dati dal programma esattamente con quel nome!!
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

{\\trasposta
tras(m)=local();
mattranspose(m);
}

\\_____________________________________________
\\ Programmi fatti per esercizio LEZ 2

\\ Lez 2 Es 0: Manipolazione di "destra". Non usa vecextract ma vector.
\\ Input: lista
\\ Output:lista con lo shift di 1 a destra
{
destra2(a)=local(n,b,ris);
n=length(a);
b=vector(n-1);
for(i=1,n-1, b[i]=a[i]);
ris=concat(a[n],b);
ris;
}

\\ Lez 2 Es 1: Sinistra
\\ Input: lista
\\ Output: shift a sinistra di un posto
{\\shift a destra di un posto
sinistra(a)=local(n,b,ris);
n=length(a);
b=vecextract(a,concat("2..",Str(n)));
ris=concat(b,a[1]);
ris;
}

\\ Lez 2 Es 1 BIS: Sinistra dalla manipolazione senza vecextract
\\ Input: lista
\\ Output: shift a sinistra di un posto
{
sinistra2(a)=local(n,b,ris);
n=length(a);
b=vector(n-1);
for(i=1,n-1, b[i]=a[i+1]);
ris=concat(b,a[1]);
ris;
}

\\ Lez 2 Es 2: Destra di n posti
\\ Input: lista e posti da shiftare
\\ Output: lista shiftata a destra di n posti
{
destran(a,posti)=local(ris);
ris=a;
for(i=1, posti,ris = destra(ris));
ris;
}

\\ Lez 2 Es 2 BIS: Sinistra di n posti
\\ Input: lista e posti da shiftare
\\ Output: lista shiftata a sinistra di n posti
{
sinistran(a,posti)=local(ris);
ris=a;
for(i=1, posti,ris = sinistra(ris));
ris;
}

\\ Lez 2 Es 3: FERMATEST
\\ Input: b,n
\\ Output: 1 se b^(n-1) = 1 mod n, 0 altrimenti
{
fermatest(b,n)=local();
if(Mod(b^(n-1),n)==1,return(1),return(0));
}

\\ Lez 2 Es 4: lista dei primi n pseudoprimi di fermat in base b
\\ Input: b, n
\\ Output: lista degli pseudoprimi di fermat sulla base b
{
fermatlista(b,n)=local();
return( vector(n,n,b^(2^(n-1))+1) );
}

\\ Lez 2 Es 5: Test Solovay-Strassen
\\sulla formula di eulero p primo allora a^((p-1)/2)) = (a/p) mod p, per ogni a.
\\Quelli che soddisfano il viceversa sono detti pseudoprimi di Eulero.
\\    b^((n-1)/2)) = (b/n) mod n     NB: n deve essere dispari!!
\\ Input: b base intera positiva, n intero positivo
\\ Output: dice se n è uno pseudoprimo di Eulero
{
eulertest(b,n)=local();
if(Mod( b^((n-1)/2 ),n)== kronecker(b,n),return(1),return(0));
}

\\ Lez 2 Es 5 BIS: LISTA PRIMI DI EULERO
\\ Input: b n
\\ Output: lista dei numeri che da 1 ad n passano il test di eulero
\\NB: i deve essere dispari!!
{
eulerlista(b,n)=local(vv);
vv=[];
for(i=1, n/2, if(eulertest(b,2*i+1)==1, vv=concat(vv,2*i+1)) );
\\Mod( b^((i-1)/2 ),i)== kronecker(b,i)
vv;
}

\\ Lez 2 Es 5 TRIS: LISTA NUMERI SI EULERO, NO PRIMI.
\\ Input: b n
\\ Output: lista dei numeri che da 1 ad n passano il test
\\         di Eulero ma non sono primi
{
eulerprimi(b,n)=local(vv,lungh,Risp);
vv=eulerlista(b,n);
lungh=length(vv);
Risp=[];
for(i=1, lungh, if(isprime(vv[i]),Risp = Risp, Risp=concat(Risp,i)   )  );
Risp;
}

\\ Lez 2 Es 6: chinese generalizzato
\\ Input: vettore di Mod(a_i,n_i) lungo k
\\ Output: Mod(a,n) che soddisfa tutti i mod del vettore
{
cinese(V) = local(Risp);
Risp=V[1];
for(i=2,length(V), Risp=chinese(Risp,V[i]));
Risp;
}

\\ Lez 2 Es 7: Generatori di (n) in Z
\\ Input: n
\\ Output: restituisce la lista dei generatori di n
{
generatori1(n) = local(Risp);
Risp=[];
for(i=1,n, if(gcd(n,i)==1, Risp=concat(Risp,i) ,Risp = Risp   ) );
Risp;
}


\\ Lez 2 Es 8: trovare i periodi di x+t in Z_17[x]/(f(x)) con f(x)
\\             irriducibile di grado 3
\\ Input: t
\\ Output: periodo di x+t
{
campofinito(t)=local(cont,ff);
cont=2;
\\irriducibile di grado 3 in Z_17:
ff =  Mod(1, 17)*x^3 + Mod(16, 17)*x^2 + Mod(15, 17)*x + Mod(6, 17);
while( Mod(x+t,ff)^cont != Mod(1,ff), cont=cont+1);
cont;
}

/*
\\ Lez 2 Es 9: MACCHINA DI TURING (facoltativo) per2
\\ Input: vettore binario del tipo [1,1,...,1] n volte.
\\ Output: riporta il vettore raddoppiato con i passaggi
\\ ad ogni step e l'elenco delle terne M(s,q)
\\ Non legge una macchina di Turing, ma è una macchina di turing

{
turing(T,nastro, q1,n1)=local(V,N,q,n);
V = T[q1,n1]; \\legge nella matrice al posto [q1,n1]
N = nastro(n1); \\legge il nastro in n1
while
 q = V[1];
 n = N(1);

}
*/

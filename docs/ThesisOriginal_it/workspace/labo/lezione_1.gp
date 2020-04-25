/*

LABORATORIO DI CALCOLO SIMBOLICO
ANNO ACCADEMICO 2010/2011

SEBASTIANO FERRARIS

\r Ferraris\esercizi

*/



\\____PRIMA LEZIONE___________________________________________________
\\____________________________________________________________________
\\____________________________________________________________________


\\_____________________________________________
\\ Programmi dati in classe LEZ 1


{\\Restituisce il primo successivo ad n (vedi foreprime)
primodopo1(n) = local(aa);
aa = n+1;
while(!isprime(aa),aa = aa+1);
aa;
}

{\\Restituisce il fattoriale di n
fattoriale(n) = local(m);
prod(m=1,n,m)
}

{\\Restituisce il fattore primo più piccolo di n (usa factor)
fattorpiccolo(n) = local();
factor(n)[1,1]
}

{\\restituisce il fattore più grande di n
fattorgrande(n) = local();
factor(n)[matsize(factor(n))[1],1]
}


\\_____________________________________________
\\ Programmi fatti per esercizio LEZ 1


\\Lez 1 Es1: Collatz semplice
\\ Input: numero iniziale di avvio calcolo sequenza
\\ Output: no. Scrive la sequenza.
{
collatz(n) = local();
while(n > 1,print(n) || if(n%2 == 0,n = n/2 , n = 3*n + 1) );
print(n);
}


\\Lez 1 ES1 bis: Collatz con la ricorsione
\\ Input: numero iniziale di avvio ricorsione
\\ Output: no. Scrive la sequenza.
{
collatzric(n) = local();
print(n);
if( n!=1, if(n%2 == 0, collatzric(n/2),collatzric(3*n + 1)), n  );
}


\\Lez 1 Es1 tris: collatz con contatore
\\ Input: numero iniziale di avvio calcolo sequenza
\\ Output: no. Scrive la sequenza con il contatore.
{
collatzcont(n) = local(cont);
cont = 0;
while(n > 1, print("count = "cont "    n = " n) ; if(n%2 == 0,n = n/2 , n = 3*n + 1) ; cont = cont + 1 )

}

\\Lez 1 Es1 ter: collatz con contatore su visualizzazione matriciale 1x2
\\ Input: numero iniziale sequenza
\\ Output: no. Stampa le matrici 1x2 con passo del calcolo e calcolo.
{
collatzmatrlong(n) = local(cont);
cont = 0;
while(n > 1,print([ cont, n ]); if(n%2 == 0, n = n/2 , n = 3*n + 1 ); cont = cont + 1   );
}

\\Lez 1 Es1 pent: collatz con contatore su matrice 1x2
\\ Input: passo richiesto del calcolo e numero iniziale sequenza
\\ Output: matrice 1x2, con passo richiesto e calcolo
{
collatzmatr(step, n) = local(cont);
cont = 0;
while(step > cont,if(n%2 == 0, n = n/2 , n = 3*n + 1 ); cont = cont + 1   );
return([ cont, n ]);
}

\\__________________________________________________________________________________
\\Lez 1 Es2: Euclide Mullin
\\ Input riceve vettore (a_i) di interi. (uso length per evitare una entrata)
\\ Output il vettore concatenato con il fattore minimo primo di Y = prod(a_i) + 1
{
euc1esercizio(V) =  local(x,y,W,minim,Risp,long);
long = length(V);
x = prod(n = 1, long, V[n]);
y = x + 1;
W = factor(y);
minim = W[1,1];
Risp = concat(V,minim);
return(Risp)
}

\\ Lez 1 Es2 bis: Euclide Mullin less local
\\ Input riceve vettore (a_i) di interi.
\\ Output il vettore concatenato con il fattore minimo primo di Y = prod(a_i) + 1
{
euc1less(V) =  local(Risp,long);
long = length(V);
Risp = concat(V,factor( prod(n = 1, long, V[n]) +1)[1,1]  ) ;
return(Risp)
}

\\ Lez 1 Es2 tris: Euclide Mullin no local
\\ Input riceve vettore (a_i) di interi.
\\ Output il vettore concatenato con il fattore minimo primo di Y = prod(a_i) + 1
{
euc1lessless(V) =  local();
return( concat(V,factor( prod(n = 1, length(V), V[n]) +1)[1,1]  )) ;
}

\\ Lez 1 Es2 ter: Euclide Mullin. Stampa dei passaggi fino a niter.
\\ Input: riceve vettore (a_i) di interi e il numero di iterazioni desiderate.
\\ Output: no. Stampa i vettori concatenati con il fattore minimo primo di Y = prod(a_i) + 1
{
euc1seq(V,niter) =  local();
print(V);
for(m = 0, niter - 1, length(V); V = concat(V,factor(prod(n = 1, length(V), V[n]) +1)[1,1]  ); print(V) );
}


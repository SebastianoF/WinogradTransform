\\----------------------------
\\ funzioni preliminari sulle matrici circolanti
\\ e funzioni ausiliarie
\\---------------------------


\\shift a destra di un posto
{
destra(a)=local(n,b,ris);
   n=length(a);
   b=vecextract(a,concat("1..",Str(n-1)));
   ris=concat(a[n],b);
ris;
}

\\shift a destra di k posti
{
destrak(a,k)=local(ris);
   ris=a;
   for(i=1,k,ris=destra(ris));
ris;
}



\\input: V vettore
\\output: V shiftato a destra di una posizione.
{
right(V)=local(Vall, Vcut, l, i);
   l = length(V);
   Vall = concat([V[l]],V);
   Vcut = vector(l, i, Vall[i]);
return(Vcut);
}

\\input: V vettore, k intero positivo
\\output: V shiftato a destra di k posizioni.
{
rightk(V,k) = local(ans, i);
   ans = V;
   for(i = 1, k, ans = right(ans));
return(ans);
}

\\input: V vettore
\\output: V shiftato a destra di una posizione.
{
left(V)=local(Vall, Vcut, l, i);
   l = length(V);
   Vall = concat(V,V[1]);
   Vcut = vector(l, i, Vall[i+1]);
return(Vcut);
}

\\input: V vettore, k intero positivo
\\output: V shiftato a destra di k posizioni.
{
leftk(V,k) = local(ans, i);
   ans = V;
   for(i = 1, k, ans = left(ans));
return(ans);
}


\\shift a sinistra di un posto
{
sinistra(a)=local();
destrak(a,length(a)-1);
}

\\riceve un array intmod e restituisce la soluzione
{
cinese(lista)=local(n,u,v);
n=length(lista);u=chinese(lista[1],lista[2]);
for(i=3,n,v=chinese(u,lista[i]);u=v);
u;
}

\\input: V vettore.
\\output: M matrice circolante generata dal vettore V
\\non usa le matrici di permutazioni.
{
circolator(V) = local(l, M, i);
   l = length(V);
   M = matrix(l,l);
   M[1,] = V;
   for(i = 2, l, M[i,] = right(M[i-1,]) );
return (M);
}

\\input: V vettore, x lemento che compare nel vettore 
\\output: vettore delle posizioni di x nel vettore V
{
position(V,x) = local(l, ans);
   l = length(V);
   ans = [];
   for(i = 1, l, 
      if(V[i] == x, 
         ans = concat(ans,i)
      )
   );
return(ans);
}


\\---------------------------
\\ permutazioni
\\---------------------------

\\composto di permutazioni
\\(2,3,1) significa 1->2, 2->3, 3->1 etc.

\\input: a,b vettori rappresentanti permutazioni.
\\output: c composizione delle permutazioni a*b
{
permcomp(a,b) = local(n,c);
   n = length(a);
   c = vector(n);
   for(i = 1, n, c[i] = a[b[i]]);
return(c);
}

\\input: n intero positivo
\\output: permutazione identica
{
permid(n) = local();
return(vector(n, j, j));
}

\\input: a permutazione
\\output: matrice corrispondente
{
perm2mat(a) = local(n, M, w, v);
n = length(a);
M = matrix(n,n);
for(j=1, n, 
   for(i=1, n, 
      if(a[i] == j, M[i,j] = 1 )
   );
);
return(M);
}


\\input: M matrice di permutazioni
\\output: premutazioni corrispondente
{
mat2perm(M) = local(l, ans);
   l = matsize(M)[1];
   ans = vector(l);
   for(i=1, l, ans[i] = position(M[i,],1)[1]);
return(ans);
}

\\input: n intero
\\output: matrice s_n
{
matshift(n)=local();
perm2mat(left(permid(n)));
}

\\input: V vettore.
\\output: M matrice circolante generata dal vettore V
\\Usa le matrici di permutazioni.
{
circolator1(V) = local(n, s, i, ans);
   n = length(V);
   s = matshift(n);
   ans = matrix(n,n);
   for(i = 0, n-1, ans = ans + V[i+1]*s^i);
return(ans);
}


\\input: a permutazione
\\output: cicli disgiunti della permutazione
{
cicli(g) = local(l, ciclo, ans, s, ini, j, cont);
l = length(g);
s = vector(l,j,j);
ans = [];
while(length(s) > 0, 
   ini = s[1];
   cont = ini;
   ciclo = [cont];
   while(g[cont] != ini,
      ciclo = concat(ciclo, g[cont]);
      cont = g[cont]
   );
   s = complemento(s,ciclo);
   ans = concat(ans,[ciclo])
   );
return(ans);
}


\\input: a permutazione
\\output: periodo di a
{
permper(a) = local(z,v,i);
   z = cicli(a);
   v = vector(length(z), i, length(z[i]));
return(lcm(v));
}

\\calcola il periodo di g (permutazione di 1..n), con cercaperiodo; si deve porre esternamente oper(a,b)=permcomp(a,b) e operid=permid(n)
{
permper2(g)=local(z,v);
cercaperiodo(g,lcm(operid));
}


\\---------------------
\\ polinomi e vettori
\\---------------------

\\inverte l'ordine in un vettore
\\usare anche Vec(Polrev(v))
{
contr(v)=local(ris,n, x);
n=length(v);ris=vector(n);forstep(x=n,1,-1,ris[n-x+1]=v[x]);
ris;
}

\\input: V vettore
\\output: V oridinato al contrario
\\usa le funzioni predefinite Vec e Polrev
{
inverter(V) = local();
return(Vec(Polrev(V)));
}


\\input: p polinomio
\\output: V vettore corrispondente a p, ordinato al contrario
\\N.B. E' l'inverso di Polrev
{
Vecrev(p)=local(ans);
   ans = inverter(Vec(p));
return(ans);
}

\\input: p polinomio
\\output: V vettore corrispondente a p, ordinato al contrario
\\N.B. Generalizza Vec
{
pol2vect(p, n) = local(m, ans);
ans = Vec(p);
m = length(ans);
if(n < m,
   return(0),
   return(inverter(concat(vector(n-m),ans)))
);
}

\\input: p1, p2, polinomi a,b in x, 
\\       di grado <= n-1, n intero positivo 
\\output: p1*p2 mod x^n-1 
{
moltpoln(a, b, n) = local(p1, p2);
   p1 = Mod(a, x^n - 1);
   p2 = Mod(b, x^n - 1);
return(lift(p1*p2));
}


\\input: M matrice circolante 
\\output: polinomio che la rappresenta
{
mat2pol(M) = local(ans);
   ans = Polrev(M[1,]);
return(ans);
}

\\input: p polinomio, l lunghezza matrice > deg(p)
\\output: matrice circolante definita dal polinomio p
{
pol2mat(p, n) = local(ans);
   ans = pol2vect(p, n);
return(circolator(ans));
}



\\---------------------
\\ matrici di Fourier
\\---------------------

\\input: n intero, p primo: n divide p-1
\\output: radice primitiva n-esima modulo p primo; 
{
primerootnp(n,p) = local(z);
   if( gcd(n,p) !=1 || !isprime(p) ,return(0));
   if( (p-1)%n != 0, return(0));
   z = znprimroot(p);
   return( z^((p-1)/n) );
}

\\matrice di Fourier nxn con n che divide p-1, p primo
{
fouriernp(n, p) = local(w, M);
   w = primerootnp(n, p);
   M = matrix(n, n, i, j, w^((i-1)*(j-1)) );
return(M);
}

\\inversa di Fourier; fournp(n,p)*Mod(circol(vettore),p)*fournpinv(n,p) è diagonale
{
fouriernpInv(n, p) = local(w, M, s);
   w = primerootnp(n,p);
   M = matrix(n, n, i, j, w^(-(i-1)*(j-1)));
   s = Mod(n,p);
reutrn(s^(-1)*M);
}


\\------------------------
\\ algebra Rn,p
\\------------------------

\\input: a,b,p polinomi in una stessa variabile, 
\\       n un intero positivo.
\\output: a*b mod p , a coeff in Zn. 
{
moltp(a,b,p,n)=local(p1,p2);
   p1 = Mod(Mod(1,n)*a,p);
   p2 = Mod(Mod(1,n)*b,p);
return( lift(lift(p1*p2)) );
}

\\come sopra in Rp,n = Zp[x]/(x^n-1)
{
moltnp(a,b,n,p) = moltp(a,b,x^n-1,p)
}

\\input: a,b polinomi in una stessa variabile, 
\\       e esponente, n un intero positivo.
\\ e può essere negativo, se a è coprimo con p
\\output: a^e mod p , a coeff in Zn. 
{
potp(a,e,p,n) = local(p1);
   p1 = Mod(Mod(1,n)*a,p);
return( lift(lift(p1^e)) );
}

\\come sopra in Rp,n = Zp[x]/(x^n-1)
{
potnp(a,e,n,p)=potp(a,e,x^n-1,p)
}

\\input: a,b polinomi in una stessa variabile, 
\\       p primo
\\output: mcd tra due polinomi a, b in Zp[x]
{
mcdp(a,b,p)=local(p1, p2);
   p1 = Mod(1,p)*a;
   p2 = Mod(1,p)*b;
return(lift(gcd(p1,p2)));
}

\\input: ap array di polinomi, mp array dei relativi moduli 
\\       p modulo. 
\\output: polinomio soluzione con il tcr
{
tcrp(ap, am, p) = local(n, u, v);
   n = length(ap);
   v = vector(n,k,Mod(Mod(1,p)*ap[k],am[k]));
   u = chinese(v);
return( lift(lift(u)) );
}

\\resto della divisione di a(x) per b(x) in Zp[x]
\\usa divrem
{
restop(a,b,p)=lift(divrem(Mod(1,p)*a,b)[2])
}

\\quoziente della divisione di a(x) per b(x) in Zp[x]
{
quozp(a,b,p)=lift(divrem(Mod(1,p)*a,b)[1])
}

{
circolmax(n,p)=local();
(lift(factormod(x^n-1,p))[,1])~
}



\\------------------------
\\ teta e inversa
\\------------------------

/*
isomorfismo tra Rp,n e il relativo prodotto di campi.
riceve un polinomio in Rp,n (ovvero un polinomio di grado < n con coeff. interi, che sono visti modulo p)
restituisce un array di s polinomi, dove s è il numero dei poli irriducibili nei quali x^n-1 si fattorizza modulo p
*/
{
teta(a,n,p)=local(f,s);
f=lift(factormod(x^n-1,p));
s=matsize(f)[1];
return(vector(s,i,restop(a,f[i,1],p)))
}

\\inverso di teta
{
tetaInv(a,n,p)=local();
tcrp(a,circolmax(n,p),p); \\teorema cinese dei resti per polinom.
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









\\----------------------------
\\ strutture algebriche del primo capitolo.
\\---------------------------


\\input: 





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
\\output: V shiftato a sinistra di una posizione.
{
left(V)=local(Vall, Vcut, l, i);
   l = length(V);
   Vall = concat(V,V[1]);
   Vcut = vector(l, i, Vall[i+1]);
return(Vcut);
}

\\input: V vettore, k intero positivo
\\output: V shiftato a sinistra di k posizioni.
{
leftk(V,k) = local(ans, i);
   ans = V;
   for(i = 1, k, ans = left(ans));
return(ans);
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

\\input: array intmod 
\\output: soluzione mediante teorema cinese dei resti
\\ generalizza chinese ad un vettore di n intmod
{
cinese(lista) = local(n, u, v);
   n = length(lista);
   u = chinese(lista[1],lista[2]);
   for(i = 3, n, 
      v = chinese(u,lista[i]); u = v
   );
return(u);
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
\\ matrici di permutazioni
\\---------------------------

\\permutazioni sono rappresentate da vettori
\\[2,3,1] significa 1->2, 2->3, 3->1 etc.

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
perm2mat(a) = local(n, M, i, j);
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

\\---------------------
\\ polinomi e vettori
\\---------------------

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
vecrev(p)=local(ans);
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


\\input: V vettore, k intero positivo
\\output: V shiftato a sinistra di k posizioni.
{
leftk(V,k) = local(ans, i);
   ans = V;
   for(i = 1, k, ans = left(ans));
return(ans);
}

\\input: V1 vettore, V2 vettore
\\output: matrice con i prodotti degli elem dei vettori mod length(V1).
{
tabProd(V1,V2) = local(r,W);
   r = length(V1);
   W = matrix(length(V2),length(V1),i,j,Mod(V1[j]*V2[i],r));
return(lift(W));
}

\\input: V1 vettore, V2 vettore
\\output: matrice con i prodotti degli elem dei vettori mod length(V1).
{
tabProd2(V1,V2,r) = local(r,W);
   r = length(V1) - 1;
   W = matrix(length(V2),length(V1),i,j,Mod(V1[j]*V2[i],r));
return(lift(W));
}

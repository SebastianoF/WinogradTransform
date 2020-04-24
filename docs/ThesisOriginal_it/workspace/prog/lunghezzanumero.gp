\\da intero a vettore dei singoli interi.

{
int2vett(n) = local(len, V);
   V = [];
   while(n != 0, V=concat(V,n%10); n=floor(n/10) );
return(length(V));
}


\\Lunghezze delle potenze di 25
\\ridà la lunghezza di 25^n
{
lunghpot25(n) = local();
return(int2vett(25^n));
}


\\ridà il vettore delle lunghezze con il numero
{
lunghpot25vet(n) = local(V);
   V=[];
   for(i=1,n, V= concat(V,lunghpot25(i) ));
return(V);
}

\\ridà il vettore delle lunghezze con il numero della posizione
{
lunghpot25vetplus(n) = local(V);
   V=[];
   for(i=1,n, V= concat(V,[[i,lunghpot25(i)]] ));
return(V);
}

\\ridà il vettore delle differenze
{
intervalvett(W) = local(V);
   V=[];
   for(i=1,length(W)-1, V= concat(V,W[i+1]-W[i] ));
return(V);
}

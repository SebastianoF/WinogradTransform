\\-------------------------------------------
\\ Funzioni di conversione
\\ in alfabeto Ascii e in alfabeto definito dall'utente.
\\-------------------------------------------

\\-----------------------------
\\ Definizione dell'alfabeto.
\\ Al primo posto (zeresimo) deve esserci lo spazio vuoto
\\ che avrà un ruolo privilegiato negli algorimi di conversione!
\\ ----------------------------


{
alfabeto = [" ","A","B","C","D","E","F","G","H","I","J","K","L",
"M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z",":", ".", "'"];
}


\\ ----------------------------
\\ Funzioni ausiliarie per conversione su
\\ alfabeto predefinito dall'utente.
\\ ----------------------------


\\input: lettera di alfabeto "A"
\\output: numero posizone della lettera partendo da zero
\\N.B. inversa di numb2letter
{
letter2numb(lettera) = local(j);
for(j=1,length(alfabeto),if(alfabeto[j]==lettera,return(j - 1)));
error("letter2numb input non valido");
}

\\input: numero corrispondente ad una lettera nell'alfabeto
\\output: lettera dell'alfabeto corrispondente dalla posizone
\\N.B. inversa di letter2numb
{
numb2letter(numero) = local();
return(alfabeto[numero + 1]);
error("numb2letter input non valido");
}

\\input: s stringa di caratteri di alfabeto
\\output: vettore di interi corrispondenti ai caratteri di s
\\come Strchr ma su alfabeto.
\\USA: Vec, letter2numb
{
string2vect(s) = local(v,i,l,aus);
   l = length(s);
   v = vector(l);
   aus = Vec(s);
   for(i = 1, l, v[i] = letter2numb(aus[i]));
return(v);
error("string2vect input non valido");
}

\\input: vettore di numeri corrispondenti
\\       a lettere di alfabeto
\\output: stringa corrispondente al vettore
\\come Vecsmall ma su alfabeto.
\\USA: numb2letter
{
vect2string(v)=local(s,l,i);
   s = "";
   l = length(v);
   for(i = 1, l, s = concat(s,numb2letter(v[i])));
return(s);
error("vect2string input non valido");
}

\\input: intero n
\\output: vettore di elementi nulli [0, ... ,0]
{
emptyVect(n) = local(V,i);
   V=vector(n);
return(V);
error("emptyVect input non valido");
}

\\input: vettore di interi V, lunghezza verme l
\\output: completa V con spazi vuoti
\\        in modo che sia di lungh multipla di l
\\USA: emptyVect
{
fullfit(V,l) = local(q,r,F);
   F = [];
   q = length(V);
   r = q%l;
   if(r > 0, F = concat(V,emptyVect(l-r)), F = V);
return(F);
error("fullfit input non valido");
}

\\input: V vettore
\\ouput: F come V privato dell'ultimo elemento
{
lastElim(V) = local(l, F, i);
   l = length(V);
   F = vector( l - 1 );
   for(i = 1, l-1, F[i] = V[i]);
return(F);
error("lastElim input non valido");
}

\\input: V vettore [1,2,3,0,0]
\\ouput: [1,2,3]
\\   elimina gli 0 corrispondenti a spazi vuoti
\\   al fondo del vettore V
\\USA: lastElim
{
trim(V) = local(l,F);
   l = length(V);
   F = V;
   while(F[l] == 0, F = lastElim(F); l = l-1);
return(F);
error("trim input non valido");
}

\\input: V vettore di interi di lunghezza multipla di l
\\output: vettore di sottovettori lunghi l, costruito da V
{
suddivector(V,l) = local(F,W,q,i,j);
   W = vector(l); F = [];
   q = length(V)/l;
   for(j=1,q,
      for(i=1,l,W[i]=V[i+(j-1)*l] );
      F = concat(F,[W]);
   );
return(F);
error("suddivector input non valido");
}

\\input: V vettore di sottovettori di pari lunghezza
\\       [[11,...,1l ],[21,...2l ], ... ,[n1,...,nl]]
\\output: vettore costruito da V come unione dei sottovettori
\\       [11, ... , nl]
{
univector(V) = local(F,l,n,i,j);
   l = length(V[1]);
   n = length(V);
   F = vector(n*l);
   for(j = 1, n,
      for(i = 1,l,F[i+(j-1)*l] = (V[j])[i] );
   );
return(F);
error("univector input non valido");
}

\\input: v singolo verme di interi in alfabeto
\\output: intero corrispondente
{
worm2int(v) = local(l,i,lAlfabeto,a);
   l = length(v);
   a = 0;
   lAlfabeto = length(alfabeto);
   for(i = 1, l, a = a*lAlfabeto + v[i]);
return(a);
error("worm2int input non valido");
}

\\input: a intero corrispondente ad un verme,
\\       l lunghezza del verme
\\output: V verme corrispondente
\\E' la pseudoinversa di worm2int
{
int2worm(a,l) = local(V, aux, lAlfabeto, i);
   V = vector(l);
   aux = a;
   lAlfabeto = length(alfabeto);
   for(i = 0, l-1,  V[l-i] = lift(Mod(aux,lAlfabeto));
                  aux = (aux - V[l-i])/lAlfabeto
   );
return(V);
error("int2worm input non valido");
}

\\input: V vettore di sottovettori di pari lunghezza
\\       [[11,...,1l],[21,...2l ], ... ,[n1,...,nl]]
\\output: [ worm2int[11,...,1l], ... , worm2int[n1,...,nl]]
\\Generalizza worm2int a tutto il vettore di vettori
{
worm2intGen(V) = local(l,i,F);
   l = length(V);
   F = vector(l);
   for(i = 1, l, F[i] = worm2int(V[i]));
return(F);
error("worm2intGen input non valido");
}

\\input: [ worm2int[11,...,1l], ... , worm2int[m1,...,ml]]
\\ vettore di interi corrispondenti a cifre codificate
\\output:V vettore di sottovettori corrispondenti
\\       [[11,...,1l],[21,...2l ], ... ,[m1,...,ml]]
\\Generalizza int2worm
{
int2wormGen(V,l) = local(m,i,F);
   m = length(V);
   F = vector(m);
   for(i = 1, m, F[i] = int2worm(V[i],l));
return(F);
error("int2wormGen input non valido");
}


\\ ----------------------------
\\ Funzioni ausiliarie per conversione alfabeto ascii
\\ che integrano le precedenti per alfabeto predefinito
\\ ----------------------------


\\input: oggetto tipo Vecsmall
\\ouput: vettore corrispondente
\\ usato nella conversione Ascii
{
vecsmall2vec(v) = local(F,i);
   F=[];
   for(i=1,length(v),F=concat(F,v[i]));
return(F);
error("Vecsmall2vec input non valido");
}

\\input: intero n
\\output: vettore di [32, 32, ... , 32]
\\        corrispondenti agli spazi vuoti " "
{
emptyVectAscii(n) = local(F,i);
   F=[];
   for(i=1,n, F = concat(F, [32]));
return(F);
error("emptyVectAscii input non valido");
}

\\input: vettore di interi V, lunghezza verme l
\\output: completa V con spazi vuoti
\\        in modo che sia di lungh multipla di l
\\USA: emptyVectAscii
{
fullfitAscii(V,l) = local(q,r,F);
   F = [];
   q = length(V);
   r = q%l;
   if(r > 0, F = concat(V, emptyVectAscii(l-r)), F = V);
return(F);
error("fullfitAscii input non valido");
}


\\input: V vettore [1,2,3,32,32]
\\ouput: [1,2,3]
\\   elimina gli 0 corrispondenti a spazi vuoti
\\   al fondo del vettore V
\\USA: lastElim
{
trimAscii(V) = local(l,F);
   l = length(V);
   F = V;
   while(F[l] == 32, F = lastElim(F); l = l-1);
return(F);
error("trimAscii input non valido");
}


\\input: v singolo verme di elementi dell'alfabeto Ascii
\\output: intero corrispondente
{
worm2intAscii(v) = local(l,i,a);
   l = length(v);
   a = 0;
   for(i = 1, l, a = a*255 + v[i]);
return(a);
error("worm2intAscii input non valido");
}

\\input: a intero corrispondente ad un verme,
\\       l lunghezza del verme
\\output: V verme corrispondente
\\E' la pseudoinversa di verme2int
{
int2wormAscii(a,l) = local(V, aux, i);
   V = vector(l);
   aux = a;
   for(i = 0, l-1,  V[l-i] = lift(Mod(aux,255));
                  aux = (aux - V[l-i])/255
   );
return(V);
error("int2worm input non valido");
}

\\input: V vettore di sottovettori di pari lunghezza
\\       [[11,...,1l],[21,...2l ], ... ,[n1,...,nl]]
\\output: [ worm2int[11,...,1l], ... , worm2int[n1,...,nl]]
\\Generalizza worm2int a tutto il vettore di vettori
{
worm2intGenAscii(V) = local(l,i,F);
   l = length(V);
   F = vector(l);
   for(i = 1, l, F[i] = worm2intAscii(V[i]));
return(F);
error("worm2intGen input non valido");
}

\\input: [ worm2int[11,...,1l], ... , worm2int[m1,...,ml]]
\\ vettore di interi corrispondenti a cifre codificate
\\output:V vettore di sottovettori corrispondenti
\\       [[11,...,1l],[21,...2l ], ... ,[m1,...,ml]]
\\Generalizza int2worm
{
int2wormGenAscii(V,l) = local(m,i,F);
   m = length(V);
   F = vector(m);
   for(i = 1, m, F[i] = int2wormAscii(V[i],l));
return(F);
error("int2wormGenAscii input non valido");
}


\\-------------------------------------------
\\ Funzioni di conversione Ascii
\\-------------------------------------------

\\input: Stringa di caratteri di Alfabeto, l lunghezza dei vermi
\\output:  vettore di interi corrispondenti ai vermi.
{
S2Vascii(S, l) = local(a, b, c, V);
   a = vecsmall2vec(Vecsmall(S));
   b = fullfitAscii(a,l);
   c = suddivector(b,l);
   V = worm2intGenAscii(c);
return(V);
error("S2Vascii input non valido");
}

\\input: V vettore di interi corrispondenti ai vermi, l lunghezza vermi
\\output: Stringa di caratteri di Alfabeto, l lunghezza dei vermi
{
V2Sascii(V, l) = local(a,b,c,S);
   a = int2wormGenAscii(V,l);
   b = univector(a);
   c = trimAscii(b);
   S = Strchr(c);
return(S);
error("V2S input non valido");
}


\\-------------------------------------------
\\ Funzioni di conversione per un alfabeto dato
\\ dall'utente
\\-------------------------------------------

\\input: Stringa di caratteri di Alfabeto, l lunghezza dei vermi
\\output:  vettore di interi corrispondenti ai vermi.
{
S2V(S, l) = local(a, b, c, V);
   a = string2vect(S);
   b = fullfit(a,l);
   c = suddivector(b,l);
   V = worm2intGen(c);
return(V);
error("S2V input non valido");
}

\\input: V vettore di interi corrispondenti ai vermi, l lunghezza vermi
\\output: Stringa di caratteri di Alfabeto, l lunghezza dei vermi
{
V2S(V, l) = local(a,b,c,S);
   a = int2wormGen(V,l);
   b = univector(a);
   c = trim(b);
   S = vect2string(c);
return(S);
error("V2S input non valido");
}


/*
Esempio di utilizzo:

? \r /.../converter.gp
? S = "MY FUNNY VALENTINE";
? S2V(S,5)
%45 = [11205201, 11740522, 1138940, 7672500]
? V2S(%,5)
%46 = "MY FUNNY VALENTINE"


? Sa = "My funny Valentine!";
? S2Vascii(Sa,5)
%48 = [327583751427, 466939396271, 411937694816, 445796842847]
? V2Sascii(%,5)
%49 = "My funny Valentine!"

*/



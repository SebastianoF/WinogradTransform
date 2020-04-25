\\ documento prof Cerruti

\\variabili globali
\\ALFABETO="abcdefghijklmnopqrstuvwxyz0123456789 ,.;:-_!?'^=()/&%$|\+*[]{}@#<>ABCDEFGHIJKLMNOPQRSTUVZ"
\\ALFABETON=Vec(Vecsmall(ALFABETO))
\\chcesare=[19, 23, 47, 20, 60, 69, 11, 83, 35, 44, 45, 85, 26, 61, 52, 39, 49, 86, 62, 70, 73, 63, 68, 7, 36, 37, 28, 71, 84, 58, 14, 64, 50, 29, 40, 10, 41, 55, 67, 9, 6, 80, 82, 66, 87, 2, 31, 1, 88, 21, 32, 8, 75, 65, 5, 74, 30, 4, 56, 22, 27, 18, 78, 24, 16, 25, 76, 34, 57, 81, 17, 48, 53, 42, 33, 54, 59, 51, 79, 72, 12, 43, 3, 46, 77, 15, 13, 38];
\\cifracesare="yT 8sz&:T{:#{T{:#{U@&:TbsTE&:zsU{A";
\\chvigenere è il testo in chiaro di cifracesare
\\cifravigenere=[59, 32, 30, 27, 20, 27, 51, 41, 57, 41, 57, 35, 41, 14, 19, 40, 30, 7, 22, 55, 23, 50, 15, 50, 3, 20, 25, 32, 37, 14, 9, 29, 9];
\\chvernam=161718;
\\cifravernam=[1098, 5426, 5863, 8128, 3654, 4739, 4531, 3074, 3659, 2136, 5476, 538, 6365, 931, 1298, 2922, 4488, 7154, 6701, 3011, 4317, 2324, 8401, 1703, 642, 6317, 5455, 6244, 2182, 219, 9934];

{\\lettera -> codice Ascii,  trasforma anche una frase: le lettere devono essere tra "".
lett2num(a)=Vec(Vecsmall(a));
}


{\\numero -> lettera; trasforma anche un vettore di numeri
num2lett(n)=Strchr(n);
}


{\\posizione della lettera a in ALFABETO
poslett(a)=posizioni(ALFABETON,lett2num(a)[1])[1];
}


{\\lettera di posizione p in ALFABETO
lettpos(p)=Strchr(ALFABETON[p]);
}


{\\atbash moderno (N.B. non si possono usare i caratteri accentati, scrivere p.e. "perche'")
atbash(testo)=local(alf,alf1,ris,n,vv);
ris=[];
alf=ALFABETO;
alfv=Vecsmall(alf);
alf1=contr(alfv);
vv=Vecsmall(testo);
n=length(vv);
for(i=1,n,g=posizioni(alfv,vv[i]);if(length(g)==0,error("il messaggio contiene un simbolo non identificato"));ris=concat(ris,alf1[g[1]]));
Strchr(ris);
}


{\\codice di cesare generalizzato: sotituzione monoalfabetica; la chiave è una permutazione dei numeri da 1 a length(ALFABETO)
cesarecod(messaggio,chiave)=local(pc,ris,alf,alfv,alfc,vv,n);
alf=ALFABETO;
vv=Vecsmall(messaggio);
n=length(vv);
pc=vecsort(chiave);
if(pc-identica(length(alf))!=vector(length(alf)),error("la chiave non è una permutazione corretta"));
ris=[];
alfv=Vecsmall(alf);
alfc=applica(chiave,alfv);
for(i=1,n,g=posizioni(alfv,vv[i]);if(length(g)==0,error("il messaggio contiene un simbolo non identificato"));ris=concat(ris,alfc[g[1]]));
Strchr(ris);
}



{\\vigenére classico; la chiave è una frase segreta;
vigenerecod(testo, chiave)=local(alf,nalf,alfv,chiavev,mm,testov,nn,ris,ch);
alf=ALFABETO;
nalf=length(alf);
alfv=Vecsmall(alf);
chiavev=Vecsmall(chiave);
mm=length(chiavev);
testov=Vecsmall(testo);
nn=length(testo);
ris=[];
ch=[];
for(i=1,mm,g=posizioni(alfv,chiavev[i]);ch=concat(ch,g[1]-1));
for(i=0,nn-1,g=posizioni(alfv,testov[i+1]);ris=concat(ris,(g[1]+ch[i%mm+1]-1)%nalf+1));
ris;
}



{\\vernam con generatore random di pari/gp; la chiave è il seme del generatore, deve essere un intero < 2^31 (in modulo)
vernamcod(testo,chiave)=local(ww,n,ris);
ww=Vecsmall(testo);
n=length(ww);
if(abs(chiave)>2^31-1,error("chiave troppo grande!"));
setrand(chiave);
ris=[];
for(i=1,n,ris=concat(ris,ww[i]+random(10^4)));
ris;
}


{\\decodifica del codice di cesare generalizzato: sotituzione monoalfabetica
cesaredecod(messaggio,chiave)=cesarecod(messaggio,inversa(chiave));
}


{\\vigenére classico; la chiave è una frase segreta;
vigeneredecod(testo, chiave)=local(nalf,alfv,chiavev,mm,testov,nn,ris,ch);
alf=ALFABETO;
nalf=length(alf);
alfv=Vecsmall(alf);
chiavev=Vecsmall(chiave);
mm=length(chiavev);
nn=length(testo);
ris=[];
ch=[];
for(i=1,mm,g=posizioni(alfv,chiavev[i]);ch=concat(ch,g[1]-1));
for(i=0,nn-1,ris=concat(ris,alfv[(testo[i+1]-ch[i%mm+1]-1)%nalf+1]));
Strchr(ris);
}


{\\vernam con generatore random di pari/gp; la chiave è il seme del generatore, deve essere in valore assoluto < 2^31
vernamdecod(testo,chiave)=local(ww,n,ris);
ww=Vecsmall(testo);
n=length(ww);
if(abs(chiave)>2^31-1,error("chiave troppo grande!"));
setrand(chiave);
ris=[];
for(i=1,n,ris=concat(ris,ww[i]-random(10^4)));
Strchr(ris);
}

{\\Applica l'iterazione di collatz per n passi
collatz(x,n)=local(a,ris);
a=x;ris=[];
for(k=1,n,if(a%2==0,a=a/2,a=3*a+1);ris=concat(ris,a));
return(ris);
}


{\\Applica l'iterazione di collatz fino a quando il risultato è 1
collatz1(x)=local(a,ris);
a=x;ris=[];
while(a!=1,if(a%2==0,a=a/2,a=3*a+1);ris=concat(ris,a));
return(ris);
}




{\\restituisce la lista delle posizioni di x nella lista lista
posizioni(lista,x)=local(nn,ris);
nn=length(lista);ris=[];
for(kk=1,nn,if(lista[kk]==x,ris=concat(ris,kk)));
ris;
}


{\\inverte l'ordine in un vettore
contr(v)=local(ris,n);
n=length(v);if(n==0,return([]));ris=vector(n);forstep(x=n,1,-1,ris[n-x+1]=v[x]);
ris;
}



{\\scambia nel vettore v gli elementi di posto i e j
scambia(v,i,j)=local(a,provv);
if(i<=0||j<=0,error("gli indici devono essere positivi!"));
if(i>length(v)||j>length(v),error("gli indici non devono essere maggiori della lunghezza del vettore!"));
a=v;
provv=a[i];
a[i]=a[j];
a[j]=provv;
a;
}


{\\rimescola un vettore v mediante m scambi
rimescola(v,m)=local(n,w);
n=length(v);
w=v;
for(i=1,m,g=random(n);h=random(n);w=scambia(w,g+1,h+1));
w;
}


{\\permutazione casuale di lunghezza n, ottenuta com rimescola 1..n, m volte
randomperm(n,m)=rimescola(vector(n,k,k),m);
}


{\\applica la permutazione p al vettore v; serve anche a comporre due permutazioni p e v
\\se la lunghezza di v è n, p deve essere una permtazione dei numeri da 1 a n
applica(p,v)=local(n,pc,w);
n=length(v);
if(length(p)!=length(v),error("le lunghezze devono essere uguali!"));
pc=vecsort(p);
if(pc-vector(n,k,k)!=vector(n),error("v non è una permutazione corretta"));
w=v;
for(i=1,n,w[i]=v[p[i]]);
w;
}


{\\permutazione inversa
inversa(p)=local(n,pc,w);
n=length(p);
pc=vecsort(p);
if(pc-vector(n,k,k)!=vector(n),error("v non è una permutazione corretta"));
w=p;
for(i=1,n,g=posizioni(p,i);w[i]=g[1]);
w;
}


{\\lista 1...n
identica(n)=vector(n,k,k);
}


{\\composto p1 o p2, si applica prima p2 e poi p1
permcomp(p1,p2)=local(ris);
nn=length(p1);
ris=vector(nn);
for(i=1,nn,ris[i]=p1[p2[i]]);
ris;
}



{\\restituisce la decomposizione di una permutazione nei suoi cicli disgiunti
permdecomp(perm)=local(nn,ris,ss,oo,kk);
nn=length(perm);
ris=[];
ss=vector(nn,k,k);
while(length(ss)!=0,ele=ss[1];oo=[ele];kk=ele;while(perm[kk]!=ele,kk=perm[kk];oo=concat(oo,kk));ris=concat(ris,[oo]);
ss=complemento(ss,oo));
ris;
}



{\\trova il periodo di una permutazione
permperiod(perm)=local(dec);
dec=permdecomp(perm);
lcm(vector(length(dec),k,length(dec[k])));
}



{\\restituisce la parità di una permutazione
permparity(perm)=local(uu,cont);
uu=permdecomp(perm);
cont=0;
nn=length(uu);
for(k=1,nn,if(length(uu[k])%2==0,cont++));
if(cont%2==1,return(1),return(0));
}



{\\crea la matrice permutazionale associata alla permutazione perm
permat(perm)=local(nn,mm);
nn=length(perm);
mm=matrix(nn,nn);
for(i=1,nn,mm[perm[i],i]=1);
mm;
}



{\\inversa di permat
matperm(a)=local(ris,n);
ris=vector(matsize(a)[1]);
n=length(ris);
for(i=1,n,ris[i]=posizioni(a[,i],1)[1]);
ris;
}



{\\media aritmetica delle componenti di un vettore
media(v)=local();
nn=length(v);
cc=sum(i=1,nn,v[i]);
cc/nn;
}



{\\dice quanti elementi y ci sono nel vettore v
quanti(v,y)=length(posizioni(v,y));
}



{\\funzione di Moebius
moeb(n)=local(ee);
if(n==1,return(1));
ee=factor(n)[,2];
if(vecmax(ee)>1,return(0));
return((-1)^length(ee));
}


{\\numirr, calcola il numero dei polinomi monici irriducibili di grado m con coefficienti in GF(q)
numirr(q,m)=local();
ris=0;
fordiv(m,d,ris=ris+moeb(m/d)*q^d);
1/m*ris;
}


{\\prodirr, calcola il prodotto dei polinomi monici irriducibili di grado m con coefficienti in GF(q)
prodirr(q,m)=local(ris);
ris=1;
fordiv(m,d,ris=ris*(x^(q^d)-x)^(moeb(m/d)));
simplify(ris);
}





{\\lista dei coniugati ad alfa di grado n su GF(q)
coniug(alfa,q,n)=vector(n,h,alfa^(q^(h-1)))
}


{\\traccia di alfa su GF(q)
traccia(alfa,q,n)=local(ee);
ee=coniug(alfa,q,n);
sum(h=1,n,ee[h]);
}


{\\norma di alfa su GF(q)
norma(alfa,q,n)=local(ee);
ee=coniug(alfa,q,n);
prod(h=1,n,ee[h]);
}


{\\Esegue due lift consecutivi
lif2(a)=lift(lift(a));
}



{\\schema dell'algoritmo potenza veloce; operid è l'elemento neutro, oper è l'operazione
pot(a,m)=local(u,t,v);
u = m; t = operid; v = a;
while(u != 0, if(u%2 == 1, t = oper(t,v)); v = oper(v,v);
u = floor(u/2));
t;
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


{\\converte n alla base b, restituisce il vettore delle cifre
converti(n,b)=local(nn,ris,gg);
if(n<0,return([]));
if(n==0,return([0]));
if(n==1,return(1));
ris=[];
nn=n;
while(nn>0, gg=divrem(nn,b); ris=concat(ris,gg[2]); nn=gg[1]);
contr(ris);
}


{\\converte n alla base b, restituisce un vettore di lunghezza L che contiene le cifra ed eventuali zeri iniziali
converti2(n,b,L)=local(nn,ris,gg);
if(n<0,return([]));
if(n==0,return(vector(L)));
if(n==1,return(concat(vector(L-1),[1])));
ris=[];
nn=n;
while(nn>0, gg=divrem(nn,b); ris=concat(ris,gg[2]); nn=gg[1]);
if(length(ris)>L, error("L troppo piccolo!"));
if(length(ris)<L,vv=vector(L-length(ris));ris=concat(ris,vv));
contr(ris);
}

{\\costruisce il campo GF(p^n)
\\utilizza un polinomio irriducibile opportuno, ff=ffinit(p,n), costruisce una radice generica al e una base
\\ {1,al,al^2,...,al^(n-1)}
\\infine restituisce l'elenco delle p^n combinazioni lineari degli elementi della base
\\con coefficienti in GF(p)
campo(p,n)=local(ff,al,base,kk);
ff=ffinit(p,n);
al=Mod(Mod(1,p)*x,ff);
base=vector(n,h,al^(h-1));
kk=[];
mm=matrix(p^n,n);
for(r=1,p^n, mm[r,]=converti2(r-1,p,n));
for(i=1,p^n,ww=Mod(1,p)*mm[i,];ww=Pol(ww);kk=concat(kk,subst(ww,x,al)));
kk;
}


{\\restituisce le n^p disposizioni con ripetizione dei numero 0..n-1, a p a p
\\in una matrice n^p x p
disprip(n,p)=local(vv,indice);
vv=matrix(n^p,p);
for(i=1,p,indice=1;for(j=1,n^(p-i),for(k=1,n,for(l=1,n^(i-1),vv[indice,i]=k-1;indice=indice+1))));
return(vv);
}



{\\vede quando il polinomio dato da ffinit è primitivo
\\è primitivo se e solo se 1 compare soltanto in una posizione (l'ultima)
\\in genere non lo è: si provi ad es. prova(2,4), prova(3,4), prova(3,3)
prova(p,n)=local(ff,al,vv,pp);
ff=ffinit(p,n);
al=Mod(x,ff);
vv=lif2(vector(p^n-1,h,al^h));
pp=posizioni(vv,1);
if(length(pp)>1,return(0),return(1));
}


{\\verifica il funzionamento di primpoly
prova2(p,n)=local(ff,al,vv);
ff=primpoly(p,n);
al=Mod(Mod(1,p)*w,ff);
vv=lif2(vector(p^n-1,h,al^h));
posizioni(vv,1);
}







{\\elimina le ripetizioni nel vettore v
elimina(v)=local();
eval(Set(v));
}


{\\restiusce il polinomio minimo in GF(p)[var] di alpha in GF(p^n)
\\il polinomio si vuole nella variabile var
\\si suppone che alpha, e quindi i suoi coniugati, siano espressi nella forma data da campo o campop
\\e che la variabile usata nei polinomi (di solito x) non si chiami XX
minimo(alpha,p,n,var)=local(cc,zar,ris);
cc=coniug(alpha,p,n);
cc=elimina(cc);
zar=Mod(1,p)*XX;
ris=lif2(prod(h=1,length(cc),zar-cc[h]));
subst(ris,XX,var);
}



{\\Restituisce in polinomio primitivo di grado n in Z_p[var]
\\L'algoritmo elimina via via gli elementi il cui ordine è un divisore proprio di p^n-1
primpoly(p, n, var) = local(s,f);
f = var^(p^n-1) -1; s = divisors(p^n -1);
for(k=1,length(s)-1, f=f/gcd(f,var^s[k] - 1));
factormod(f, p)[1, 1];
}



{\\costruisce il campo GF(p^n)
\\utilizza un polinomio primitivo opportuno, ff=primpoly(p,n,var), nella variabile desiderata var,  e costruisce una radice generica al
\\restituisce l'insieme delle P^n-1 potenze distinte di al
\\unito con 0, al primo posto. All'ultimo posto c'è 1
campop(p,n,var)=local(ff,al,kk);
ff=primpoly(p,n,var);
al=Mod(Mod(1,p)*var,ff);
kk=vector(p^n-1,h,al^h);
kk=concat([0],kk);
}


{
\\riceve tre polinomi a,b,p in una stessa variabile e calcola a*b mod (p,n) , dove n è un intero positivo
moltp(a,b,p,n)=local(y,z);
y=Mod(Mod(1,n)*a,p);z=Mod(Mod(1,n)*b,p);
return(lift(lift(y*z)))}


{
\\riceve due polinomi in una stessa variabile e due interi e=esponente, n=modulo
\\n deve essere positivo, e può essere negativo, se a è coprimo con p
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



{\\matrice delle tracce Tr(a_i*a_j), dove lista={a_1,...,a_n}
\\lista è una base se e solo se il determinante di mattracce è diverso da 0
\\si noti che gli elementi della matrice sono nel campo base
mattracce(lista,p,n)=local(mm);
mm=matrix(n,n);
for(i=1,n,for(j=1,n,mm[i,j]=traccia(lista[i]*lista[j],p,n)));
mm;
}



{\\matrice le cui colonne sono i coniugati di a_1,...,a_n
matpowers(lista,p,n)=local(mm);
mm=matrix(n,n);
for(i=1,n,for(j=1,n,mm[i,j]=lista[j]^(p^(i-1))));
mm;
}


{\\trova il minimo t tale che p^t>=b, se b>1. Se b==1, restituisce 0
trovat(b,p)=local(t);
if(b==1,return(0));
t=1;
while(p^t<b,t++);
t;
}


{\\calcola polcyclo(n) con la fuzione di Moebius e i polinomi x^d-1
polciclo(n)=local(ris);
ris=1;
fordiv(n,d,ris=simplify(ris*(x^d-1)^moeb(n/d)));
ris;
}


{\\cerca l'ordine di un polinomio mod p
cercaordine(pol,p)=local(tt,ss,nn,vv);
tt=factormod(pol,p);
ss=matsize(tt);
nn=ss[1];
vv=vector(nn);
for(i=1,nn,xx=Mod(x,tt[i,1]);id=xx^0;
fordiv(p^poldegree(tt[i,1])-1,n,if(xx^n==id,vv[i]=n)));
for(i=1,nn,vv[i]=vv[i]*(p^trovat(tt[i,2],p)));
lcm(vv);
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


{\\calcola i laterali ciclotomici, dati q ed n
laterali(n,q)=local(g);
g=generato(n,q);
classi(n,g);
}


{\\calcola il termine n-esimo della sequenza ricorrente che ha polinomio caratteristico f(x) e dati iniziali nel vettore v
ricor(f,v,n)=local(mm,vv,mat);
mm=length(v);
if(mm!=poldegree(f),error("il polinomio ha grado errato!"));
if(n<=mm-1,return(v[n+1]));
mat=matcompanion(f);
vv=v*mat^(n-mm+1);
vv[mm];
}


{\\calcola (modulo m) il termine n-esimo della sequenza ricorrente che ha polinomio caratteristico f(x) e dati iniziali nel vettore v
ricormod(f,v,m,n)=local(mm,vv,mat);
mm=length(v);
if(mm!=poldegree(f),error("il polinomio ha grado errato!"));
v=Mod(v,m);
if(n<=mm-1,return(v[n+1]));
mat=matcompanion(f);
mat=Mod(mat,m);
vv=v*mat^(n-mm+1);
lift(vv[mm]);
}


{\\test di Fermat
fermatest(b,n)=local(g);
g=Mod(b,n);
if(lift(g^(n-1))==1,1,0);
}


{\\lista degli pseudoprimi di Fermat
fermatlista(b,ini,fine)=local(g);
if(ini%2==0,ini=ini+1);
g=[];
forstep(k=ini,fine,2,if(fermatest(b,k)==1 && isprime(k)==0,g=concat(g,k)));
g;
}


{\\test di miller
millertest(a,n)=local(aa,bb,ss,tt);
if(n%2==0,return(0));
if(gcd(a,n)!=1,return(0));
ss=0;tt=n-1;
while(tt%2==0,tt=tt/2;ss++);
aa=Mod(a,n);
bb=aa^tt;
if(bb==Mod(1,n),return(1));
for(k=0,ss-1,if(bb==Mod(-1,n),return(1));bb=bb^2);
return(0);
}


{\\lista degli pseudoprimi di Miller
millerlista(b,ini,fine)=local(g);
if(ini%2==0,ini=ini+1);
g=[];
forstep(k=ini,fine,2,if(millertest(b,k)==1 && isprime(k)==0,g=concat(g,k)));
g;
}

{\\test di Eulero
eulertest(b,n)=local(g,h);
if(gcd(b,n)!=1,return(0));
g=Mod(b,n);
h=Mod(kronecker(b,n),n);
if(g^((n-1)/2)==h,1,0);
}



{\\lista degli pseudoprimi di Eulero
eulerlista(b,ini,fine)=local(g);
if(ini%2==0,ini=ini+1);
g=[];
forstep(k=ini,fine,2,if(eulertest(b,k)==1 && isprime(k)==0,g=concat(g,k)));
g;
}



{\\test AKS
akstest(a,r,n)=local(xx,bb);
if(gcd(n,r)!=1,return(0));
if(gcd(n,a)!=1,return(0));
xx=Mod(Mod(1,n)*x,x^r-1);
bb=Mod(a,n);
return((xx+bb)^n==xx^n+bb);
}


{\\lista degli pseudoprimi AKS
akslista(a,r,ini,fine)=local(g);
if(ini%2==0,ini=ini+1);
g=[];
forstep(k=ini,fine,2,if(akstest(a,r,k) && !isprime(k),g=concat(g,k)));
g;
}

{\\verifica, data una lista di interi n che sono aks-pp rispetto all'esponente r (base qualsiasi), che r divide sempre n^2-1
aksver(r,lista)=local(nn);
nn=length(lista);
for(i=1,nn,if((lista[i]^2-1)%r!=0,return(i)));
1;
}


{\\calcola ricormod(x^2-h*x+k,[a,b],m,n)
ricmod2(a,b,h,k,m,n)=local(z,y);
z=Mod([0,1;-k,h],m);
y=z^n;
return(lift(a*y[1,1]+b*y[1,2]))}


{\\test per gli pp di fibonacci di prima specie
fib1test(n)=local();if(n%5==0,return(0));
if(ricmod2(0,1,1,-1,n,n)==kronecker(5,n)%n,return(1),return(0))}



{\\test per gli pp di fibonacci di seconda specie
fib2test(n)=local();if(n%5==0,return(0));
if(ricmod2(0,1,1,-1,n,n-kronecker(5,n))==0,return(1),return(0))}



{\\lista degli pp di Fibonacci di prima specie
fib1lista(ini,fine)=local(g);
if(ini%2==0,ini=ini+1);
g=[];
forstep(k=ini,fine,2,if(fib1test(k) && !isprime(k),g=concat(g,k)));
g;
}


{\\lista degli pp di Fibonacci di seconda specie
fib2lista(ini,fine)=local(g);
if(ini%2==0,ini=ini+1);
g=[];
forstep(k=ini,fine,2,if(fib2test(k) && !isprime(k),g=concat(g,k)));
g;
}


{\\intersezione delle due liste precedeni
fiblista(ini,fine)=local(ris);
ris=[];
if(ini%2==0,ini=ini+1);
forstep(n=ini,fine,2,if(fib1test(n)==1 && fib2test(n)==1 && !isprime(n),ris=concat(ris,n)));
return(ris)}


{\\applica fib1test e fib2test insieme
fibtest(n)=local(zz,mm,kk);
if(n%5==0, return(0));
zz=Mod([0,1;1,1],n);
mm=zz^n;
kk=kronecker(5,n);
if(kk==1,if(mm[1,1]==Mod(0,n) && mm[1,2]==Mod(kk,n),return(1),return(0)));
if(mm[2,2]==Mod(0,n) && mm[1,2]==Mod(kk,n),return(1),return(0));
}



{\\fibtest generalizzato
fibgentest(P,Q,n)=local(xx,jj,cc,delta);
if(n==2,return(1));
if(gcd(n,2*Q)!=1,return(0));
xx=Mod(1,n)*[0,1;-Q,P];
delta=P^2-4*Q;
if(n<=delta && isprime(n),return(1));
jj=kronecker(delta,n);
if(jj==0,return(0));
if(jj==1,cc=lift(xx^(n-1));if(cc[1,2]%n==0 && (cc[2,2]-1)%n==0,return(1),return(0)));
cc=lift(xx^n);if(cc[2,2]%n==0 && (cc[1,2]+1)%n==0,return(1),return(0));
}



{\\come fiblista, ma più veloce
fiblistav(ini,fine)=local(g);
if(ini%2==0,ini=ini+1);
g=[];
forstep(k=ini,fine,2,if(fibtest(k)==1 && isprime(k)==0,g=concat(g,k)));
g;
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


{\\matrice circolante -> polinomio che la rappresenta
mat2pol(a)=local();
Polrev(a[1,]);
}

{\\inversa di mat2pol
pol2mat(b,n)=local(c);
c=pol2vetc(b,n);
circol(c);
}



{\\dice se la matrice cicolante nxn a è invertibile in Q[x]/(x^n-1)
matinvQ(a)=local(n);
n=matsize(a)[1];
if(gcd(Polrev(a[1,]),x^n-1)==1,1,0);
}


{\\fattorizza f(x) in GF(p^e)=[x]
\\es: factoringfq(x^17-1,2,4)
factoringfq(a,p,e)=local(f);
f=lift(ffinit(p,e,t));
lift(lift(factorff(a,p,f)));
}


{\\trova i gradi dei polinomi irriducibili
\\nel cui prodotto si fattorizza x^n-1 in GF(q)[x]
gradi(n,q)=local(z);
z=classi(n,generato(n,q));
s=length(z);
v=vector(s,j,length(z[j]));
vecsort(v);
}


{\\restituisce l'ordine del gruppo moltiplicativo di GF(q)[x]/(x^n-1)
circolgrorder(n,q)=local(v);
v=gradi(n,q);
prod(i=1,length(v),q^(v[i])-1);
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
\\riceve tre polinomi a,b,p in una stessa variabile e calcola a*b mod (p,n) , dove n è un intero positivo
moltp(a,b,p,n)=local(y,z);
y=Mod(Mod(1,n)*a,p);z=Mod(Mod(1,n)*b,p);
return(lift(lift(y*z)))}

{
\\riceve due polinomi in una stessa variabile e due interi e=esponente, n=modulo
\\n deve essere positivo, e può essere negativo, se a è coprimo con p
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
\\somma di polinomi nell'algebra Rp,n = Zp[x]/(x^n-1)
sumnp(a,b,n,p)=local(aa,bb);
aa=Mod(Mod(1,p)*a,x^n-1);bb=Mod(Mod(1,p)*b,x^n-1);
lift(lift(aa+bb));
}


{
\\potenza di un polinomio nell'algebra Ap,n
potnp(a,e,n,p)=potp(a,e,x^n-1,p)}


{
\\tcrp : riceve un array di polinomi, un array dei relativi moduli e il modulo p e restituisce il polinomio soluzione
tcrp(ap,am,p)=local(n,u,v);
n=length(ap);v=vector(n,k,Mod(Mod(1,p)*ap[k],am[k]));u=chinese(v);
return(lift(lift(u)))}



{
\\quoziente della divisione di a(x) per b(x) in Zp[x]
quozp(a,b,p)=lift(divrem(Mod(1,p)*a,b)[1])}

{
\\isomorfismo tra Rp,n e il relativo prodotto di campi.
\\riceve un polinomio in Ap,n (ovvero un polinomio di grado < n con coeff. interi, che sono visti modulo p)
\\restituisce un array di s polinomi, dove s è il numero dei poli irriducibili nei quali x^n-1 si fattorizza modulo p
teta(a,n,p)=local(f,s);
  f=lift(factormod(x^n-1,p));
  s=matsize(f)[1];
return(vector(s,i,restop(a,f[i,1],p)))
}

{
\\inverso di teta
gamm(a,n,p)=local();
tcrp(a,circolmax(n,p),p);
}

\\in entrambe l'ordine della fattorizzazione non è quello dato dal vettore L ma
\\ è quello dato dalla fattorizzazione implementata in pari gp.

\\trasforma un polinomio in un elemento dell'algebra R(n,p)
polpol(pol,n,p)=Mod(Mod(1,p)*pol,x^n-1)


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

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\/////////////////////////////////////////////////////////////////

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

{
\\vettore casuale di lunghezza n con componenti comprese tra 0 e m-1
randvett(n,m)=vector(n,k,random(m));
}

{
\\probabilità che una circolante nxn abbia il determinante divisibile per il primo p, n generico
npprobgen(n,p)=local();
if(n%p==0,return(0));
1.-circolgrorder(n,p)/(p^n);
}

{
\\probabilità che una circolante nxn abbia il determinante divisibile per il primo p, posto che n divida p-1
npprob(n,p)=local();
if((p-1)%n != 0,return(0));
1.-(1-1/p)^n;
}

{
\\probabilità ricavata sperimentalmente
npprobex(n,p,casi)=local();
favorevoli=0;
for(i=1,casi,if(matdet(circol(randvett(n,p)))%p==0,favorevoli++));
1.*favorevoli/casi;
}

{
\\log in base 2
log2(a)=log(a)/log(2);
}

{
\\riceve n
\\restituisce un fattore di n
\\o un r coprimo con r tale che il periodo moltiplicativo di n mod r sia > [4log(n)^2]
akspasso2(n)=local(r,L,v);
r=2;L=floor(4*log2(n)^2);
while(r<n,d=gcd(r,n);if(d>1,return(d));v=znorder(Mod(n,r));
if(v>L,break());
r=r+1);
[r,v,L,isprime(r),n,floor(16*log(n)^5)]
}

{\\AKS
aks(n)= local(w,r,L,l,d,v,xx,bb);
if(n==2,return("primo"));
\\termina se n è pari diverso da 2
if(n%2==0,print("n è pari");return(0));
\\verifica se è una potenza perfetta
w=ispower(n);
if(w>0,print("n è una potenza perfetta");return(w));
\\trova r
r=2;L=floor(4*(log2(n))^2);
while(r<n,d=gcd(r,n);if(d!=1,print("n è composto");return(d));v=znorder(Mod(n,r));
if(v>L,break());
r=r+1);
if(r==n, return("primo"));
\\loop finale
l=2*(floor(sqrt(r)))*floor(log2(n))+1;
xx=Mod(x,x^r-1);
for(b=1,l,
bb=Mod(b,n);
if((xx+bb)^n!=xx^n+bb,print("n è composto");return([r,bb])));
return(1);
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

\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
\\restituisce la matrice T che rappresenta tetaQ
del(n)=local(ff,poli,s,T,resti);
ff=factor(x^n-1);
poli=ff[,1]~;
s=length(poli);
T=matrix(n,n);resti=vector(n);
for(j=1,n,resti[j]=[];for(i=1,s,resti[j]=concat(resti[j],pol2vetc(divrem(x^(j-1),poli[i])[2],poldegree(poli[i])))));
for(i=1,n,T[,i]=mattranspose(resti[i]));
T;
}

{
\\applica del(n) come pura trasformazione lineare
dellin(a)=local(n);
n=length(a);
mattranspose(del(n)*mattranspose(a));
}

{
\\applica del(n) come biiezione tra Q[x]/(x^1-1) e il prodotto dei Q[x]/(polciclotomico fattore di x^n-1)
\\è la formattazione corretta di dellin(a,n)
\\dà lo stesso risultato di tetaQ(a,n)
delpol(a,n)=local(ap,b,ff,poli,s,v,d);
ap=pol2vetc(a,n);
b=dellin(ap);
ff=factor(x^n-1);
poli=ff[,1]~;
s=length(poli);
v=vector(s);ini=0;
for(k=1,s,d=poldegree(poli[k]);v[k]=Polrev(vector(d,j,b[ini+j]));ini+=d);
v;
}

{
\\inversa di del(n)
gamlin(a)=local(n,mm);
n=length(a);
mm=del(n)^-1;
mattranspose(mm*mattranspose(a));
}

{
\\dà lo stesso risultato di del(n)^-1, ma costruisce la matrice senza bisogno di conoscere del(n)
gam(n)=local(ff,poli,s,ris,ma,d,ini);
ff=factor(x^n-1);
poli=ff[,1]~;
s=length(poli);
ris=matrix(n,n);ini=1;
ma=matid(s);
for(k=1,s,d=poldegree(poli[k]);ris[ini,]=pol2vetc(tcrQ(ma[k,],poli),n);for(j=ini+1,ini+d-1,ris[j,]=destra(ris[j-1,]));ini+=d);
mattranspose(ris);
}

{
\\restituisce la matrice T che rappresenta teta
delnp(n,p)=local(ff,poli,s,T,resti);
if(gcd(n,p)!=1,return(0));
ff=lift(factormod(x^n-1,p));
poli=ff[,1]~;
s=length(poli);
T=matrix(n,n);resti=vector(n);
for(j=1,n,resti[j]=[];for(i=1,s,resti[j]=concat(resti[j],pol2vetc(restop(x^(j-1),poli[i],p),poldegree(poli[i])))));
for(i=1,n,T[,i]=mattranspose(resti[i]));
T;
}

{
\\dà lo stesso risultato di delnp(n,p)^-1, ma costruisce la matrice senza bisogno di conoscere delnp(n,p)
gamnp(n,p)=local(ff,poli,s,ris,ma,d,ini);
ff=lift(factormod(x^n-1,p));
poli=ff[,1]~;
s=length(poli);
ris=matrix(n,n);ini=1;
ma=matid(s);
for(k=1,s,d=poldegree(poli[k]);ris[ini,]=pol2vetc(tcrp(ma[k,],poli,p),n);for(j=ini+1,ini+d-1,ris[j,]=destra(ris[j-1,]));ini+=d);
mattranspose(ris);
}


{\\diagonalizza a blocchi una circolante a nxn con i coefficienti mod p
circdiag(a,n,p)=lift(Mod(1,p)*delnp(n,p)*mattranspose(a)*gamnp(n,p));
}


{\\diagonalizza a blocchi una circolante a nxn con i coefficienti in Q
circdiagQ(a,n)=del(n)*mattranspose(a)*gam(n);
}


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\/////////////////////////////



{
provaecm(nn,ee)=local(xx,yy,aa,bb,elli);
xx=random(nn);yy=random(nn);aa=random(nn);
bb=yy^2-xx^3-aa*xx;
elli=ellinit([0,0,0,Mod(aa,nn),Mod(bb,nn)]);
ellpow(elli,[xx,yy],ee);
}


{
hadamard(n)=local(hh,mh,m1,m2);
hh=[1,1;1,-1];
if(n==1,return(hh));mh=hh;
for(j=1,n-1,m1=concat(mh,mh);m2=concat(mh,-mh);mh=concat(mattranspose(m1),mattranspose(m2)));
mh;
}


{
	ellido(nn)=local(xx,yy,aa,bb);
	xx=random(nn);yy=random(nn);aa=random(nn);
	bb=yy^2-xx^3-aa*xx;
	[ellinit([0,0,0,Mod(aa,nn),Mod(bb,nn)]),[xx,yy]];
}


{
ellisum(a,P,Q)=local(lam, mu, x3, x1, y1, x2, y2);
 if(P == [0], return(Q)); if(Q == [0], return(P)); x1 = P[1]; x2 = Q[1]; y1 = P[2];
 y2 = Q[2];if(x1==x2&&y1==-y2,return([0]));
 if(P == Q, lam = (3*x1^2 + 2*a[2]*x1 + a[4] - a[1]* y1)/(2* y1 + a[1]* x1 + a[3]), lam = (y2 - y1)/(x2 - x1));
 mu = y1 - lam *x1; x3 = lam^2 + a[1]* lam - a[2] - x1 - x2;
 [x3, -(lam + a[1])* x3 - mu - a[3]];
}


{
ellipot(a,P,m)=local(u,t,v);
u = m; t = [0]; v = P;
while(u != 0, if(u%2 == 1, t = ellisum(a,t,v)); v = ellisum(a,v,v);
u = floor(u/2));
t;
}


{
ellisumn(z,P,Q,n)=local(lam, mu, x3, x1, y1, x2, y2,den,fat,a);
 if(P == [0], return(Q)); if(Q == [0], return(P));
 a=z*Mod(1,n);
 x1 = Mod(P[1],n); x2 = Mod(Q[1],n); y1 = Mod(P[2],n);
 y2 = Mod(Q[2],n);if(x1==x2&&y1==-y2,return([0]));
 if(P == Q, den=(2* y1 + a[1]* x1 + a[3]); fat=gcd(n,lift(den)); if(fat>1 && fat!=n,return(fat));
 lam = (3*x1^2 + 2*a[2]*x1 + a[4] - a[1]* y1)*den^-1,den= (x2 - x1); fat=gcd(n,lift(den)); if(fat>1 && fat!=n,return(fat));
 lam = (y2 - y1)*den^-1);
 mu = y1 - lam *x1; x3 = lam^2 + a[1]* lam - a[2] - x1 - x2;
 lift([x3, -(lam + a[1])* x3 - mu - a[3]]);
}


{
ellipotn(a,P,m,n)=local(u,t,v);
u = m; t = [0]; v = P;
while(u != 0, if(u%2 == 1, t = ellisumn(a,t,v,n);if(type(t)=="t_INT",return(t))); v = ellisumn(a,v,v,n);if(type(v)=="t_INT",return(v));
u = floor(u/2));
lift(t);
}


{
trovaexp(p,B)=local(ee);
if(p>B,return(0));
ee=1;
while(p^ee<=B,ee++);
ee-1;
}


{
ecm(n,B)=local(xx,yy,cc,ex,aa,bb,pun,ris);
cc=1;
while(cc==1,
xx=random(n);yy=random(n);aa=random(n);
bb=(yy^2-xx^3-aa*xx)%n;
gg=gcd(4*aa^3+27*bb^2,n);
if(gg>1 && gg!=n, return(gg));
if(gg==1,pun=[xx,yy];cc=0));
forprime(pr=2,B,ex=trovaexp(pr,B);
for(j=1,ex,ris=ellipotn([0,0,0,aa,bb],pun,pr,n);if(ris==[0],return([0]));if(type(ris)=="t_INT",return(ris));pun=ris));
error("fallimento!");
}



{
	ellpunti1(a,b,p)=local(e,i,j,punti);
e=ellinit([0,0,0,a,b]);
punti=[];
for(i=0, p-1, for(j=0, p-1, q=[Mod(i,p), Mod(j,p)];
                            if(ellisoncurve(e,q)==1, punti=concat(punti, [q]))));
return(lift(punti))}



{
	ellpunti2(a,b,p)=local(e,i,x,s,punti,g,q1,q2);
e=ellinit([0,0,0,a,b]);
punti=[];
for(i=0, p-1, x=Mod(i,p);
              s=x^3+a*x+b;
              if(kronecker(lift(s),p)==1,g=s^(1/2); q1=[x,g];
                                          q2=[x,-g];
                                          punti=concat(punti,[q1]);
                                          punti=concat(punti,[q2]));
              if(lift(s)==0, q=[x, s]; punti=concat(punti,[q])));
return(lift(punti))}

{
	ellord(z,p)=p+1-ellap(ellinit(z),p)
}

\\ [0, 0, 0, -191, -510] [4] r=0
\\ [1, -1, 1, -122, 1721] [12] r=0
\\  [1, 0, 0, -49423080, 130545230400] [2,8] r=2

msp=785963102379428822376694789446897396207498568951;
msa=317689081251325503476317476413827693272746955927;
msb=79052896607878758718120572025718535432100651934;
msgx=771507216262649826170648268565579889907769254176;
msgy=390157510246556628525279459266514995562533196655;

{convlist(a)=local(n,p,q,pa,pb,qa,qb,lp,lq);n=length(a);pb=1;pa=a[1];lp=[pa];qb=0;qa=1;lq=[qa];for(i=2,n,p=a[i]*pa+pb;q=a[i]*qa+qb;lp=concat(lp,p);
lq=concat(lq,q);pb=pa;pa=p;qb=qa;qa=q);return([lp,lq])}

{conv(a,k)=local(b);b=convlist(contfrac(a,k));return(b)}

{conv1(a,k)=local(c,n,lis);lis=[];c=conv(a,k);n=length(c[1]);for(i=1,n,lis=concat(lis,c[1][i]/c[2][i]));return(lis)}

{\\risolve, se possibile, x^2 = a (mod n)
qsolven(a,p)=local(unmezzo,h,k);
if(gcd(a,p)!=1,return(gcd(a,p)));
if(kronecker(a,p)!=1,return(0));
if(p==2,return(1));
unmezzo=(p+1)/2;
h=1;k=a;
while(kronecker(h^2-4*k,p)!=-1,h=h+1);
ris=ricmod(2,h,h,k,(p+1)/2,p);
ris=(unmezzo*ris)%p;
ris=min(ris,p-ris);
if(ris^2%p==a,return(ris),return(0));
}


{fconv(a)=local(n,p,q,pa,pb,qa,qb,lp,lq);n=length(a);pb=1;pa=a[1];lp=[pa];qb=0;qa=1;lq=[qa];for(i=2,n,p=a[i]*pa+pb;q=a[i]*qa+qb;lp=concat(lp,p);
lq=concat(lq,q);pb=pa;pa=p;qb=qa;qa=q);return([lp,lq])}

{
prodd(d,a,b)=[a[1]*b[1]+d*a[2]*b[2],a[1]*b[2]+a[2]*b[1]]
}

{rad(d) = local(a, P, Q, s, v);
    s = sqrt(d); v = floor(s);
    a = []; P = [0]; Q = [1];
    k = 1;
    while(1, a = concat(a, floor((P[k] + v)/Q[k]));
      P = concat(P, a[k]*Q[k] - P[k]);
      Q = concat(Q, (d - P[k + 1]^2)/Q[k]);
      if(a[k] == 2*v, break); k=k+1);
    return(a)}

    {pell(D) = local(r, n, a, v, b);
    r = rad(D); n = length(r) - 1;a = fconv(r);
    if(n%2 == 0,
      return([a[1][n], a[2][n]]),
      return([[a[1][n], a[2][n]], prodd(D,[a[1][n], a[2][n]],[a[1][n], a[2][n]])]));
     }

     {radtutti(d) = local(a, P, Q, s, v);
    s = sqrt(d); v = floor(s);
    a = []; P = [0]; Q = [1];
    k = 1;
    while(1, a = concat(a, floor((P[k] + v)/Q[k]));
      P = concat(P, a[k]*Q[k] - P[k]);
      Q = concat(Q, (d - P[k + 1]^2)/Q[k]);
      if(a[k] == 2*v, break); k=k+1);
    return([a,P,Q])}


{provarad(d,n) = local(a, P, Q, s, v);
    s = sqrt(d); v = floor(s);
    a = []; P = [0]; Q = [1];
    k = 1;
    for(g=1,n, a = concat(a, floor((P[k] + v)/Q[k]));
      P = concat(P, a[k]*Q[k] - P[k]);
      Q = concat(Q, (d - P[k + 1]^2)/Q[k]);
      if(a[k] == 2*v, break); k=k+1);
    return(a)}



{pell2(D) = local(r, n, a, v, b);
    r = rad(D); n = length(r) - 1;
    if(n%2 == 0, a = fconv(r);
      return([a[1][n], a[2][n]]), a = fconv(r);
      v = concat(r, vecextract(r, "2..")); b = fconv(v);
      return([[a[1][n], a[2][n]], [b[1][2*n], b[2][2*n]]]))}


{
      maplus(k)=[floor(sqrt(k)),k;1,floor(sqrt(k))]
}


{
      maminus(k)=[floor(sqrt(k)+1),k;1,floor(sqrt(k)+1)]
}

{
gauss(a)=1/(a-floor(a))
}

{
gaussm(a)=1/(ceil(a)-a)
}

{
iterplus(k,n)=local(ris,vv);
ris=matrix(2,n);
vv=[1,0]~;for(i=1,n,vv=maplus(k)*vv;ris[1,i]=vv[1];ris[2,i]=vv[2]);
ris;
}


{
iterminus(k,n)=local(ris,vv);
ris=matrix(2,n);
vv=[1,0]~;for(i=1,n,vv=maminus(k)*vv;ris[1,i]=vv[1];ris[2,i]=vv[2]);
ris;
}

{fconvm(a)=local(n,p,q,pa,pb,qa,qb,lp,lq);n=length(a);pb=1;pa=a[1];lp=[pa];qb=0;qa=1;lq=[qa];for(i=2,n,p=a[i]*pa-pb;q=a[i]*qa-qb;lp=concat(lp,p);
lq=concat(lq,q);pb=pa;pa=p;qb=qa;qa=q);return([lp,lq])}

{per2(a,b)=2*a/b+a^2}

{per2m(a,b)=-2*a/b+a^2}

{
contfracm(a,n)=local(it,b,ris);
it=[a];b=a;for(i=1,n,b=gaussm(b);it=concat(it,b));
ris=[ceil(a)];for(i=1,n,ris=concat(ris,ceil(it[i+1])));
ris;
}

{
charplus(k)=2*floor(sqrt(k))/(k-floor(sqrt(k))^2)
}

{
gaussp(a,p)=local();
(a-truncate(a+O(p)))^-1;
}

{
padicint(a,p)=local();
truncate(a+O(p));
}

\\{
\\padicserlist(a,p)=local();


\\test basato sulla ERH
{
erhtest(n)=local(B);
B=2*log(n)^2;
forprime(pr=2,B,if(gcd(n,pr)>1,return(0));if(eulertest(pr,n)==0,return(0)));
return(1);
}

\\\\
\\\\\

{\\generalizzazione della formula di Kesava Menon
kesava(n,q)=local(gg,zz);
if(q==0,gg=coprimi(n);zz=length(gg);return(1/zz*sum(h=1,zz,gcd(gg[h]-1,n))));
gg=generato(n,q);
zz=length(gg);
1/zz*sum(h=1,zz,gcd(gg[h]-1,n))
}




{
stestm(n)=local(al,be,a1,a2);
if(n==2||n==3||n==5||n==7||n==11||n==13,return(1));
if(gcd(n,2310)!=1,return(0));
if(n%13==0,return(0));
al=Mod(Mod(1,n)*x,x^3-x^2-4*x-1);
be=al^n;
a1=x^2-2*x-2;a2=-x^2+x+3;
if(be-x==0||be-a1==0||be-a2==0,return(1));
return(0);
}




{\\shanks primality test, vedere le mie ricerche
stestc(n)=local(rr,al,be,a1,a2);
if(n==2||n==3||n==5||n==7||n==11||n==13,return(1));
if(gcd(n,2310)!=1,return(0));
rr=n%13;
if(rr==0,return(0));
al=Mod(Mod(1,n)*x,x^3-x^2-4*x-1);
be=al^n;
a1=x^2+(n-2)*x-2;a2=(n-1)*x^2+x+3;
if((rr==1||rr==5||rr==8||rr==12),if((be-x==0),return(1),return(0)));
if((rr==2||rr==3||rr==10||rr==11),if((be-a1==0),return(1),return(0)));
if((rr==4||rr==6||rr==7||rr==9),if((be-a2==0),return(1),return(0)));
}




{
stest(n)=local(al);
al=Mod(Mod(1,n)*x,x^3-x^2-4*x-1);
lift(lift(al^n));
}





{\\Dati p primo > 2 ed e > 0 calcola i p^e-1 quadrati ortogonali risultanti
\\dal campo finito
fflatin(p,e)=local(ff,mm,ris);
ff=campop(p,e,x);
mm=matrix(p^e,p^e);
ris=vector(p^e-1);
for(f=2,p^e,forvec(x=[[1,p^e],[1,p^e]],mm[x[1],x[2]]=posizioni(ff,ff[f]*ff[x[1]]+ff[x[2]])[1]);ris[f-1]=mm);
ris;
}



{retta(m,b,x)=m*x+b}



{\\Dati p primo > 2 ed e > 0 calcola i p^e-1 quadrati ortogonali risultanti
\\dalla geometria affine generata dal campo finito
fglatin(p,e)=local(mm,ff,ris);
mm=matrix(p^e,p^e);
ff=campop(p,e,x);
ris=vector(p^e-1);
for(m=2,p^e,for(b=1,p^e,forvec(x=[[1,p^e],[1,p^e]],if(retta(ff[m],ff[b],ff[x[1]])==ff[x[2]],mm[x[1],x[2]]=b)));ris[m-1]=mm);
ris;
}


{\\test di primalità basato su Frobenius
frtest(n)=local();
y=Mod(Mod(1,n)*x,x^2-x-1);
z=lift(lift(y^(2*n)-y^n-1));
if(z==0,return(1));
0;
}


{\\prodotto di Kronecker di matrici
kron(a,b)=local(la,lb,lc,c);
la=length(a);lb=length(b);lc=la*lb;c=matrix(lc,lc);
forvec(x=[[0,lc-1],[0,lc-1]],c[x[1]+1,x[2]+1]=a[floor(x[1]/lb)+1,floor(x[2]/lb)+1]*b[x[1]%lb+1,x[2]%lb+1]);
c;
}


{
millerfiblista(a,b)=local(ris);
ris=[];
forstep(h=a,b,2,if(millertest(2,h)==1&&fibgentest(2,2,h)==1&&!isprime(h),
ris=concat(ris,h)));
ris;
}



{
genfiblista(a,b)=local(ris);
ris=[];
forstep(h=a,b,2,if(fibgentest(2,2,h)==1&&!isprime(h),
ris=concat(ris,h)));
ris;
}


{
isover(n)=local(f);
f=factor(n)[,1];
h=length(f);
a=znorder(Mod(2,n));
for(i=1,h,if(znorder(Mod(2,f[i]))!=a,return(0)));
1;
}


{
overum(n)=local(k);
k=znorder(Mod(2,n));
for(i=1,k-1,if(gcd(n,2^i-1)!=1,return(0)));
1;
}


{
overumlista(a,b)=local(ris);
ris=[];
forstep(h=a,b,2,if(overum(h)==1&&!isprime(h),
ris=concat(ris,h)));
ris;
}


{\\elenco delle etichette
etic(n,p)=local(v,d,et);
if(n%p==0,return(0));
v=laterali(n,p);d=length(v);et=[];
for(i=1,d,et=concat(et,vecmin(v[i])));
vecsort(et);
}


\\{\\matrice delta compressa
\\le(n,p)=local(et,m,pol,xx,be,z);
\\et=etic(n,p);z=length(et);m=znorder(Mod(p,n));pol=primpoly(p,m,x);xx=Mod(x,pol);be=xx^((p^m-1)/n);
\\matrix(z,z,i,j,be^(et[i]*et[j]));
\\}


{\\group - inverse nel prodotto di campi
kinv(a,n,p)=local(zz,ris,m);
zz=circolmax(n,p);m=length(zz);
ris=vector(m);
for(i=1,m,if(a[i]!=0,ris[i]=lift(lift((Mod(Mod(1,p)*a[i],zz[i]))^-1))));
ris;
}


{\\group - inverse nel prodotto di campi, caso Q
kinvQ(a,n)=local(zz,ris,m);
zz=factor(x^n-1);m=matsize(zz)[1];
ris=vector(m);
for(i=1,m,if(a[i]!=0,ris[i]=lift((Mod(a[i],zz[i,1]))^-1)));
ris;
}


{\\group-inverse in R(n,q)
ginv(a,n,p)=gamm(kinv(teta(a,n,p),n,p),n,p)
}



{\\group-inverse su Q
ginvQ(a,n)=gammaQ(kinvQ(tetaQ(a,n),n),n)
}



{\\group-inverse di una cicolante, caso n,p
circinv(a,n,p)=circol(pol2vetc(gamm(kinv(teta(Polrev(a[1,]),n,p),n,p),n,p),n));
}


{\\group-inverse di una cicolante, caso Q
circinvQ(a,n)=circol(pol2vetc(gammaQ(kinvQ(tetaQ(Polrev(a[1,]),n),n),n),n));
}


{\\polinomio al contrario
polcontr(a,n,p)=moltnp(Pol(pol2vetc(a,n)),1,n,p);
}


{\\polinomio coniugato
polinv(a,n,p)=moltnp(x,polcontr(a,n,p),n,p);
}


{\\polinomio reciproco
polrec(a,n,p)=moltnp(x^poldegree(a,x),polinv(a,n,p),n,p);
}


{\\restituisce l'idempotente e che genera l'ideale (a)
polidem(a,n,p)=local(g,h,d,e,e2,lam);
g=mcdp(a,x^n-1,p); if(poldegree(g,x)==0,return(1));
h=quozp(x^n-1,g,p);
d=bezout(Mod(1,p)*g,Mod(1,p)*h);
d=lift(d);
e=g*d[1];
e2=moltnp(e,e,n,p);
lam=quozp(e2,e,p);
quozp(e,lam,p);
}

{\\peso di un polinomio
polw(a)=local(h,m,t);
if(a==0,return(0));
h=Vec(a);m=length(a); t=posizioni(h,0);
m-length(t);
}

{\\distanza tra due polinomi
pold(a,b)=polw(a-b);
}

{\\distanza tra due polinomi in R(n,p)
poldnp(a,b,n,p)=polw(sumnp(a,-b,n,p));
}


{\\distanza minima in un sottoinsieme di R(n,p)
distnp(s,n,p)=local(m,d);
m=length(s);d=n;
forvec(t=[[1,m],[1,m]],d=min(d,poldnp(s[t[1]],s[t[2]],n,p)),2);
d;
}

{\\dati a in R(n,p) e b in R(n,q) applica il tcr agli array dei coefficienti e ricostruisce il polinomio
fromnpq(a,b,n,p,q)=local(za,zb,vv);
za=pol2vetc(a,n);zb=pol2vetc(b,n);
vv=lift(chinese(Mod(za,p),Mod(zb,q)));
Polrev(vv);
}

{\\dati a in R(n,p) e b in R(n,q) applica il tcr agli array dei coefficienti e ricostruisce il vettore
fromnpqv(a,b,n,p,q)=local(za,zb,vv);
za=pol2vetc(a,n);zb=pol2vetc(b,n);
vv=lift(chinese(Mod(za,p),Mod(zb,q)));
vv;
}


{\\inversa di fromnpq
tonpq(y,n,p,q)=local(zz,z1,z2);
zz=pol2vetc(y,n);z1=zz%p;z2=zz%q;
[Polrev(z1),Polrev(z2)];
}


{\\usa n=7,p=11,q=13
cod1(a,b)=local(m,chiave,c11,cc12,c21,c22,z1,z2);
m=143;

chiave=[x+x^2,4*x+1,11+6*x];

c11=gamm([a,1,1],7,11);
c12=gamm([1,b,chiave[1]],7,11);

c21=gamm([1,1,chiave[2],1],7,13);
c22=gamm([1,1,1,chiave[3]],7,13);

z1=fromnpqv(c11,c21,7,11,13);
z2=fromnpqv(c12,c22,7,11,13);

[z1,z2];


}




{\\usa n=7,p=11,q=13
cod12(a,b)=local(m,chiave,c11,cc12,c21,c22,z1,z2);
m=143;

chiave=[x+x^2,4*x+1,11+6*x];

c11=gamm([a,1,1],7,11);
c12=gamm([1,b,chiave[1]],7,11);

c21=gamm([1,1,chiave[2],1],7,13);
c22=gamm([1,1,1,chiave[3]],7,13);

z1=fromnpqv(c11,c21,7,11,13);
z2=fromnpqv(c12,c22,7,11,13);

[z1,z2];

}

{\\p deve essere congruo a 4 mod 7, a,b,c sono visti mod p
cod2(a,b,c,p)=gamm([random(p),a+b*x+c*x^2,random(p)],7,p);
}

{\\diagonalizza c+c~ su Q; c=circol(a)
symdiag(a)=circdiagQ(circol(a)+mattranspose(circol(a)),length(a));
}


{\\produce una matrice circolante simmetrica con m valori iniziali ini
\\se z=0 la matrice è di ordine 2m
\\se z=1 la matrice è di ordine 2m-1
symcircol(a,z)=local(m,b);
m=length(a);
if(z==0,return(circol(concat(a,contr(a)))),b=vecextract(a,vector(m-1,h,h));return(circol(concat(vecextract(a,vector(m,h,h)),contr(b)))));
}

{\\periodic autocorrelation function of a family with m rows, of lengh n
corr(g)=local();
m=matsize(g)[1];
n=matsize(g)[2];
ris=vector(n-1);
for(s=1,n-1,ris[s]=sum(i=0,n-1,sum(k=1,m,g[k,i+1]*g[k,(i+s)%n+1])));
ris;
}


{\\autocorrelation
autocorr(a)=local(u);
u=matrix(2,length(a));
u[1,]=a;u[2,]=a;
corr(u);
}


{
pp2=[1,0,0,0,0,0;0,0,0,0,0,1;0,0,0,0,1,0;0,0,0,1,0,0;0,0,1,0,0,0;0,1,0,0,0,0];
}































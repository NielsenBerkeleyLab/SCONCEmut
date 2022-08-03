#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> /*remember to use -lm flag when compiling*/
#include <time.h>
# include "ranlib.h"

//-------------DATA STRUCTURES AND GLOBALS-------------
double GENOMELENGTH; //This is made as a global for easy access in all functions.  It should not change doing the execution of the program
int root,tip,numleafs,comma=0; /*globals used to read in the tree*/
int SNPNUMBER; //counts the total number of SNPs
#define Cmax 1000 //This is the maximum number of chromosomes allowed
FILE *infile;

/*This is used for making linked lists of chromosomes to store a whole genome*/
struct Genome {
  struct Genome* nextChr; /*points to next chromosome*/
  struct Chr* chr; /*pointer to first node in the chromsome it points to*/
};

struct Genome** AllGenomes;

/*This is used for making linked lists of events on chromosomes to store a whole chromosome*/
struct Chr {
  double pos; /*position of the event*/
  int type;  /*type of event. 0: start, 1: end, 2: SNV*/
  struct Chr* next;
};

struct allele_count {
  int cell;
  int derived;
  int ancestral;
  struct allele_count* next;
};

/*This is used to count SNPs*/
struct SNP {
  double pos; /*position of SNP*/
  int snpnumber;
  struct SNP* next;
};

/*this is used for storing the tree*/
struct treenode {
  int child[2];
  int parent;
  double bl;
};
struct treenode *tree;

struct param {
  int numbins;
  int negbin_r;
  double number_reads;
  double readlength;
  double errorrate;
  double *binrates;
};




//----------------------------------------------

/*ERRROR MESSAGE FUNCTION*/
void myerror(int i)
{
  printf("Internal error number %i\nTerminating program\n",i);
  exit(-1);
}

/*minimum and maximum functions*/
double max(double a, double b)
{
  if (a>b) return a;
  else return b;
}

double min(double a, double b)
{
  if (a<b) return a;
  else return b;
}

//------------------------------FUNCTIONS NEEDED  FOR SETTING UP THE TREE AND READING THE TREEFILE--------

/*subfunction needed by ‘getclade*/
void linknodes(int i,int j,int nodee) /*linking i parent to nodee and j parent to nodee*/
{
  tree[nodee].child[0]=j;
  tree[nodee].child[1]=i;
  tree[i].parent=nodee;
  tree[j].parent=nodee;
}

/*subfunction needed by ‘getclade*/
int specsearch()
{
  char ch;
  int i=1;
  ch = fgetc(infile);
  if ((ch!=')')&&(ch!='(')&&(ch!=',')&&(ch!=' ')&&(ch!='\t')&&(ch!='\n')&&(ch!=EOF)){
    ungetc(ch, infile);
    fscanf(infile,"%d",&tip);
    while ((ch=fgetc(infile))!=':'&&(i<10))
      i++;
    tree[tip+numleafs-2].child[0]=-1;
    tree[tip+numleafs-2].child[1]=-1;
    while ((ch=(fgetc(infile)))==' ');
    ungetc(ch,infile);
    fscanf(infile,"%lf",&tree[tip+numleafs-2].bl);
    //printf("\nbranchlength of node %i =%f",tip+numleafs-1,tree[tip+numleafs-2].bl);
    return 1;
  }
  else {
    ungetc(ch, infile);
    return 0;
  }
}

/*subfunction needed by ‘getclade*/
int getnodenumb()
{
  char c;
  int i,j=0;
  fpos_t position;
  i=0;
  fgetpos(infile, &position);
  do{
    c=fgetc(infile);
    if (c==',')
      i++;
    if (c=='(')
      j=j-1;
    if (c==')')
      j++;
  }
  while ((j<0)&&(c!=EOF));
  fsetpos(infile,&position);
  return (i+comma+1);
}

/*some old code for reading a Newick tree*/
int getclade()
{
  int n1, n2, n3;
  char ch;

  do{
    if (specsearch()==1){
      /*tip++;*/
      return tip+(numleafs-1);
    }
    ch = fgetc(infile);
    if (ch==','){
      comma++;
    }
    if (ch==')'){
      if ((ch=fgetc(infile))!=':'){
        ungetc(ch,infile);
      }
      else {
        do{
          ch=(fgetc(infile));
        }while((ch=='\n')||(ch==' '));
        ungetc(ch,infile);
        fscanf(infile,"%lf",&tree[n3-1].bl);
      }
      return n3;
    }
    if (ch=='('){
      n3=getnodenumb();
      n1=getclade();
      n2=getclade();
      linknodes(n1-1,n2-1,n3-1);
    }
  }
  while (ch!=';');
  return -1;
}

//------------------------------------------------------------------------------


//--------------FUNCTIONS NEEDED FOR RV GEBNEATIONS---------------------



int static z_rndu=137;
void SetSeed(int seed)
{
  z_rndu = 170*(seed%178) + 137;
}


double uniform()
{
  /*

     U(0,1): AS 183: Appl. Stat. 31:188-190
     Wichmann BA & Hill ID.  1982.  An efficient and portable
     pseudo-random number generator.  Appl. Stat. 31:188-190

     x, y, z are any numbers in the range 1-30000.  Integer operation up
     to 30323 required.

     Suggested to me by Z. Yang in the Lord's year of 1995 who also provided me with
     the source code used here.
     */
  static int x_rndu=11, y_rndu=23;
  double r;

  x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
  y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
  z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
  if (x_rndu<0) x_rndu+=30269;
  if (y_rndu<0) y_rndu+=30307;
  if (z_rndu<0) z_rndu+=30323;
  r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
  return (r-(int)r);
}

/* Generating an exponential*/

double expo(double c)

{
  return - (1.0/c)*log(uniform());
}



/*  LOCALGAMMLN - Natural log of the gamma function.    */

double localgammln(xx)
  double xx;
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

#define PI 3.141592653589793

double binomialdev(pp,n)
  double pp;
  int n;
{
  double localgammln();
  double uniform();
  int j;
  static int nold=(-1);
  double am,em,g,angle,p,bnl,sq,t,y;
  static double pold=(-1.0),pc,plog,pclog,en,oldg;

  p=(pp <= 0.5 ? pp : 1.0-pp);
  am=n*p;
  if (n < 25) {
    bnl=0.0;
    for (j=1;j<=n;j++)
      if (uniform() < p) ++bnl;
  } else if (am < 1.0) {
    g=exp(-am);
    t=1.0;
    for (j=0;j<=n;j++) {
      t *= uniform();
      if (t < g) break;
    }
    bnl=(j <= n ? j : n);
  } else {
    if (n != nold) {
      en=n;
      oldg=localgammln(en+1.0);
      nold=n;
    } if (p != pold) {
      pc=1.0-p;
      plog=log(p);
      pclog=log(pc);
      pold=p;
    }
    sq=sqrt(2.0*am*pc);
    do {
      do {
        angle=PI*uniform();
        y=tan(angle);
        em=sq*y+am;
      } while (em < 0.0 || em >= (en+1.0));
      em=floor(em);
      t=1.2*sq*(1.0+y*y)*exp(oldg-localgammln(em+1.0)
          -localgammln(en-em+1.0)+em*plog+(en-em)*pclog);
    } while (uniform() > t);
    bnl=em;
  }
  if (p != pp) bnl=n-bnl;
  return bnl;
}

#undef PI


//------------------------------------------------------------------------------

//--------------------UTILITIES NEEDED FOR MANIPULATING THE GENOME IN VARIOUS WAYS------------------



/*This is used to insert a new node in a linked list defining a chromosome*/
void insertNodeAfterChr(struct Chr* prev_node, double pos, int type)
{
  if (prev_node == NULL) myerror(3);
  struct Chr* new_node =(struct Chr*) malloc(sizeof(struct Chr));
  new_node->pos = pos;
  new_node->type = type;
  new_node->next = prev_node->next;
  prev_node->next = new_node;
}

/*This is used to insert a new node in a linked list defining a genome*/
void insertNodeAfterGenome(struct Genome* prev_node, struct Chr* firstchr)
{
  if (prev_node == NULL)  myerror(4);
  struct Genome* new_node =(struct Genome*) malloc(sizeof(struct Genome));
  new_node->chr = firstchr;
  new_node->nextChr = prev_node->nextChr;
  prev_node->nextChr = new_node;
}

/*This is used to insert a new node in the first position in linked list defining a genome*/
struct Genome* insertNodeBeforeGenome(struct Genome* next_node, struct Chr* firstchr)
{
  if (next_node == NULL)  myerror(5);
  struct Genome* new_node =(struct Genome*) malloc(sizeof(struct Genome));
  new_node->nextChr = next_node;
  new_node->chr = firstchr;
  return new_node;
}

//CALCULATES THE LENGTH OF A CHROMOSOME
double lengt_of_chromosome(struct Chr* n)
{
  double p, L =0.0;

  while (n != NULL) {
    if (n->type==0) p=n->pos;
    else if (n->type==1) L += n->pos-p;
    n = n->next;
  }
  return L;
}

//PRINTS A CHROMOSOME
void printChromosome(struct Chr* n)
{
  double L;
  L = lengt_of_chromosome(n);
  while (n != NULL) {
    printf(" (%.1lf, %d) ", n->pos,n->type);
    n = n->next;
  }
  printf(" Length %lf\n",L);
}

//PRINTS A WHOLE GENOME
void printGenome(struct Genome* n)
{
  while (n != NULL) {
    printChromosome(n->chr);
    n = n->nextChr;
  }
  printf("\n");
}

/*makes a new chromosome of length 'length'*/
/*chromosome length is measured in doubles from 0 to 'length'*/
void makeCHR_no_gaps(struct Genome* G, int length)
{
  struct Chr* start  = NULL;
  struct Chr* stop = NULL;

  start = (struct Chr*)malloc(sizeof(struct Chr));
  stop = (struct Chr*)malloc(sizeof(struct Chr));

  G->chr = start;
  start->pos = 0; // assign data in first node
  start->type = 0;
  start->next = stop; // Link first node with second

  stop->pos = length; // assign data to second node
  stop->type = 1;
  stop->next = NULL;
}

/*make a new diploid genome of length 'length'*/
void makeDiploidGenomeNoGaps(struct Genome* first, double length)
{
  struct Genome* second = NULL;

  second = (struct Genome*)malloc(sizeof(struct Genome));

  makeCHR_no_gaps(first, length);
  first->nextChr = second; // Link first node with second

  makeCHR_no_gaps(second, length);
  second->nextChr = NULL;
}

/*make a new duplicated diploid genome of length 'length'*/
void makeDiploidGenomeNoGapsDuplicated(struct Genome* first, double length)
{
  struct Genome* second = NULL;
  struct Genome* third = NULL;
  struct Genome* fourth = NULL;

  second = (struct Genome*)malloc(sizeof(struct Genome));
  third = (struct Genome*)malloc(sizeof(struct Genome));
  fourth = (struct Genome*)malloc(sizeof(struct Genome));

  makeCHR_no_gaps(first, length);
  first->nextChr = second; // Link first node with second

  makeCHR_no_gaps(second, length);
  second->nextChr = third; // Link first node with third

  makeCHR_no_gaps(third, length);
  third->nextChr = fourth; // Link first node with fourth

  makeCHR_no_gaps(fourth, length);
  fourth->nextChr = NULL;
}


/*copies the chromosome 'n' and returns pointer to the new copied chromosome*/
struct Chr* copyChr(struct Chr* n)
{
  if (n == NULL) return NULL;
  else {
    struct Chr* newNode  = (struct Chr*)malloc(sizeof(struct Chr));
    newNode->pos = n->pos;
    newNode->type = n->type;
    newNode->next = copyChr(n->next);
    return newNode;
  }
}

/*copies the genome 'g' and returns pointer to the new copied genome*/
struct Genome* copyG(struct Genome* g)
{
  if (g == NULL) return NULL;
  else {
    struct Genome* newNode  = (struct Genome*)malloc(sizeof(struct Genome));
    newNode->chr = copyChr(g->chr);
    newNode->nextChr = copyG(g->nextChr);
    return newNode;
  }
}

/*this function is used exclusively by make_amplification*/
/*it copies nodes and makes a new terminal node*/
struct Chr* copyChr_with_max(struct Chr* n, double end)
{
  if (n == NULL) return NULL;
  struct Chr* newNode  = (struct Chr*)malloc(sizeof(struct Chr));
  if (n->pos > end)
  {
    newNode->pos = end;
    newNode->type = 1;
    newNode->next = NULL;
  }
  else {
    newNode->pos = n->pos;
    newNode->type = n->type;
    newNode->next = copyChr_with_max(n->next, end);
  }
  return newNode;
}

/*make amplification of chromosome n from post 'start' to position 'end' to form a new chromosome*/
/*this function allows 'end' to be larger than the length of the chromosome but start must be shorter than the length*/
struct Chr* make_amplification(struct Chr* n, double start, double end)
{

  if (n == NULL) myerror(41);
  else if (start < n->pos) return copyChr_with_max(n, end);
  else {
    while (n != NULL && n->pos <= start)
      n = n->next;
    if  (n== NULL) myerror(2);
  }
  struct Chr* startNode  = (struct Chr*)malloc(sizeof(struct Chr));
  startNode->pos = start;
  startNode->type = 0;
  startNode->next = copyChr_with_max(n, end);
  return startNode;
}

void make_SNV(struct Chr* n, double pos)
{

  if (n == NULL) myerror(-1);
  else if (n->pos > pos) myerror(-10);
  while (n->next != NULL && n->next->pos < pos)  {
    n = n->next;
    if  (n->next== NULL) myerror(-11);
  }
  struct Chr* newNode  = (struct Chr*)malloc(sizeof(struct Chr));
  newNode->pos = pos;
  newNode->type = 2;
  newNode->next = n->next;
  n->next = newNode;
}


/*finds and removes an empty chromosome*/
struct Genome* clean_genome(struct Genome* G)
{
  struct Genome* temp = G;
  struct Genome* prev = NULL;//have to do this because the list is not doubly linked
  struct Genome* node;

  while (G != NULL)
  {
    if (G->chr == NULL)
    {
      node = G;
      if (prev==NULL) temp = G->nextChr;
      else prev->nextChr = G->nextChr;
      G = G->nextChr;
      free(node);
    }
    else
    {
      prev = G;
      G = G->nextChr;
    }
  }
  return temp;
}

/*this deletes all nodes in chr after the initial node called*/
void delete_chr(struct Chr* c)
{
  if (c->next != NULL) delete_chr(c->next);
  free(c);
}

/*this deletes all nodes in genome after the initial node called*/
void delete_genome(struct Genome* G)
{
  delete_chr(G->chr);
  if (G->nextChr != NULL) delete_genome(G->nextChr);
  free(G);
}



/* used exclusively by make_deletion to delete nodes and set up the new starting node where the deletion ends*/
struct Chr* deleteChr_with_(struct Chr* n, double end)
{
  struct Chr* node;
  while (n != NULL && n->pos <= end) {
    node = n;
    n = n->next;
    free(node);
  }
  if (n==NULL || n->type==0) return n;
  else {
    struct Chr* newNode  = (struct Chr*)malloc(sizeof(struct Chr));
    newNode->pos = end;
    newNode->type = 0;
    newNode->next = n;
    return newNode;
  }
}


/*make deletion of chromosome n from post 'start' to position 'end'*/
/* the start position must be less than the end of the chromosome*/
struct Chr* make_deletion(struct Chr* n, double start, double end)
{

  //printf("makes deletion between %lf and %lf\n",start,end);

  struct Chr* node=n;
  /*first checking that we have a chromosome which should consist of at least two node.  These first two lines can probably safely be deleted if there are no other errors in the program*/
  if (n == NULL) myerror(6);
  if (n->next == NULL) myerror(7);
  /*If deletion starts before the chromosome start.  May be allowable under some models*/
  if (n->pos >= start) return deleteChr_with_(n, end);
  else if (n->next->pos <= start) {
    while (n->next != NULL && n->next->pos <= start)
      n = n->next;
    if  (n->next== NULL) myerror(8);
  }
  /*if deletion starts in deleted segment.  May be allowable under some models*/
  if (n->type == 1){
    n->next = deleteChr_with_(n, end);
    return node;/*this indicates that the first node is still the same*/
  }
  /*if deletion starts in an actual existing segment of DNA*/
  else {
    struct Chr* startNode  = (struct Chr*)malloc(sizeof(struct Chr));
    startNode->pos = start;
    startNode->type = 1;
    startNode->next = deleteChr_with_(n->next, end);
    n->next=startNode;
    return node;/*this indicates that the first node is still the same*/
  }

}

int number_of_chrmosomes(struct Genome* n)
{
  int i=0;

  while (n != NULL) {
    i++;
    n = n->nextChr;
  }
  return i;
}

/*assumes an array of  dubles of length equal to number of chromosomes has been allocated*/
double store_genome_length(struct Genome* n, double *GL)
{
  int i=0;
  double L=0;

  while (n != NULL) {
    GL[i]=lengt_of_chromosome(n->chr);
    L += GL[i];
    i++;
    n = n->nextChr;
  }
  return L;
}

//function for testing if the genome is internally consistent
void test_chromosome(struct Chr* n)
{
  int state;

  if (n == NULL) myerror(100);
  state = 0;
  do {
    if (n->type == 0)
    {
      if (state==1) myerror(101);
      else state = 1;
    }
    else if (n->type == 1)
    {
      if (state==0) myerror(102);
      else state = 0;
    }
    else if (n->type ==2)
    {
      if (state==0) myerror(103);
    }
    else myerror(104);
    if (n->next==NULL && n->type != 1) myerror(105);
    n = n->next;
  } while (n != NULL);
}

void test_genome(struct Genome* n)
{

  while (n != NULL) {
    test_chromosome(n->chr);
    n = n->nextChr;
  }
}

//finds a position in a genome relative to reference genome
double findpoint(struct Chr* n, double point)
{
  double lastpos, sum = 0.0;

  if (n==NULL) myerror(22);
  if (point==0.0) return 0.0;
  lastpos = n->pos;
  n = n->next;
  while (n != NULL) {
    if (n->type!=0) sum += n->pos - lastpos;
    if (sum >= point) return n->pos - sum + point;
    else if (n->type!=1) lastpos = n->pos;
    n = n->next;
  }
  return GENOMELENGTH;
}
//------------------------------------------------------------------------------


//-----------------CORE FUNCTIONS SIMULATING GENOME EVOLUTION ALONG A LINEAGE-------------------


struct Genome* mutate_genome(struct Genome* G, double *L, double r_sum, double startrate,  double rate, double r_del, double r_amp, double l_del, double l_amp, double theta)
{
  double startpoint, cumu, u, U;
  struct Genome* temp = G;
  int i;

  U = rate*uniform();
  cumu = 0.0;
  i=0;
  while (U > cumu)
  {
    cumu += startrate;
    if (U<cumu) /*the event happens at the start of the chromosome*/
    {/*a deletion*/
      if ((u=uniform())<r_del/r_sum)
      {
        G->chr = make_deletion(G->chr, 0.0, findpoint(G->chr,expo(1.0/l_del)));
        if (G->chr==NULL) temp = clean_genome(temp);
      }
      else {/*an amplification*/
        temp = insertNodeBeforeGenome(temp, NULL);
        temp->chr = make_amplification(G->chr, 0.0, findpoint(G->chr,expo(1.0/l_amp)));
      }
    }
    else {
      cumu += r_sum*L[i];
      if (U<cumu)
      {/*the event happens interior in the chromosome*/
        startpoint = L[i]*uniform();
        if ((u=uniform())<r_del/r_sum)
        {/*a deletion*/
          G->chr = make_deletion(G->chr, findpoint(G->chr,startpoint), findpoint(G->chr,startpoint+expo(1.0/l_del)));
          if (G->chr==NULL) temp = clean_genome(temp);
        }
        else if (u<(r_del+r_amp)/r_sum){/*an amplification*/
          temp = insertNodeBeforeGenome(temp, NULL);
          temp->chr = make_amplification(G->chr, findpoint(G->chr,startpoint), findpoint(G->chr,startpoint+expo(1.0/l_amp)));
        }
        else {
          make_SNV(G->chr, findpoint(G->chr,startpoint));}
      }
    }
    i++;
    G = G->nextChr;
  }
  //  printf("After mutation we obtain the following genome\n");
  //  printGenome(temp);
  return(temp);
}



/*this function simulates evolution along an edge in the phylogeny*/
/*
r_del: rate of deletion per unit of genome
r_amp: rate of amplification per unit of genome
l_del: mean length of deletions
l_amp: mean length of amplifications
T: branch length
so for a genome of length 200 the instantaneous rate of deletion is 200*r_del
*/

struct Genome* simulate_edge(struct Genome* G, double r_del, double r_amp, double l_del, double l_amp, double theta, double T)
{
  int k;
  double L[Cmax], TL, rate, startrate, r_sum, t = 0.0;

  do {
    k = number_of_chrmosomes(G); /*this part can be made much faster using a linked list - no need to recalculate everything*/
    if (k>=Cmax) myerror(1);
    TL = store_genome_length(G, L); /*this part can be made much faster using a linked list - no need to recalculate everything*/
    r_sum = r_del + r_amp + theta;
    startrate = r_del*l_del+r_amp*l_amp;
    rate = r_sum*TL + (double)k*startrate;//this is the trick to ensure a uniform rate across the genome
    t += expo(rate);
    //printf("Mutation at time %lf, rate: %lf, # of chromosomes: %i, until T: %lf, TL: %lf\n",t,rate,k, T,TL);
    if (t<T) G = mutate_genome(G, L,  r_sum,  startrate,   rate, r_del,  r_amp,  l_del,  l_amp, theta);
    //printGenome(G);
    //test_genome(G);
  } while (t<T);
  return G;
}

//------------------------Functions for printing out bed files and simulating conditional on the true states in each cell-----------------------------


double count_in_bin(struct Genome* G, struct param parameters, double bins[])
{
  double total, binsize, prev_pos, end, add;
  int i, nb;
  struct Chr* c;

  nb = parameters.numbins;
  binsize = GENOMELENGTH/(double)nb;
  total=0.0;
  while (G != NULL)
  {
    c = G->chr;
    while (c != NULL)
    {
      if (c->type!=0)
      {
        end = c->pos;
        for (i=0; i<nb; i++)
        {
          add =  min(end, (double)(i+1)*binsize) - max(prev_pos, (double)i*binsize);
          //printf("bin: (%.3lf, %.3lf), interval: (%.3lf, %.3lf): %.3lf\n",(double)i*binsize,(double)(i+1)*binsize, prev_pos, end, add);
          if (add > 0){
            bins[i] += add;
            total += add;
          }
        }
      }
      prev_pos=c->pos;
      c = c->next;
    }
    G = G->nextChr;
  }
  return total;
}

struct SNP*  screen_print_SNPs(struct SNP* mySNPS)
{
  struct SNP* temp = mySNPS;

  printf("SNPS:\n");
  while (mySNPS != NULL){
    if (mySNPS->pos > 0.0)
      printf("%i: %lf, ",mySNPS->snpnumber, mySNPS->pos);
    mySNPS = mySNPS->next;
  }
  printf("\n");
  return temp;
}

struct SNP* find_SNPS(int cellnumber, struct Genome* G, struct SNP* mySNPS)
{
  struct Chr* c;
  struct SNP* temp = mySNPS;

  //  printf("Cell number %i, ",cellnumber);
  while (G != NULL)
  {
    c = G->chr;
    while (c != NULL && c->type!= 1)
    {
      if (c==NULL) myerror(-17);
      else if (c->type==2) {
        mySNPS=temp;
        while (mySNPS->next != NULL && mySNPS->next->pos<=c->pos)
          mySNPS=mySNPS->next;
        if (mySNPS->pos==c->pos);
        else
        {
          struct SNP* newnode  = (struct SNP*)malloc(sizeof(struct SNP));
          SNPNUMBER++;
          newnode->pos = c->pos;
          newnode->next = mySNPS->next;
          mySNPS->next=newnode;
          // printf("insert SNP at pos %lf after pos %lf (%lf)\n",c->pos,mySNPS->pos,newnode->pos);
        }
        //screen_print_SNPs(mySNPS);
      }
      c = c->next;
    }
    G = G->nextChr;
  }
  return temp;
}


struct SNP* enumrate_SNPS(struct SNP* mySNPS)
{

  struct SNP* temp = mySNPS;
  int i=0;

  while (mySNPS != NULL)
  {
    mySNPS->snpnumber=i-1;
    mySNPS=mySNPS->next;
    i++;
  }
  return temp;
}

struct SNP* get_allelel_counts_in_bins(int cellnumber, struct Genome* G, struct SNP* mySNPS, int numbins, int ***allelecount)

{

  int i, step, bin = 0;
  struct Chr* c;
  struct SNP* temp = mySNPS;
  double start;

  while (G != NULL)
  {
    c = G->chr;
    start=c->pos;
    mySNPS=temp->next;
    while (c != NULL)
    {
      if (c->type==2) {
        while (mySNPS != NULL && mySNPS->pos<=c->pos)
        {
          if (mySNPS->pos>start && mySNPS->pos != c->pos) allelecount[cellnumber][mySNPS->snpnumber][0]++;
          else if (mySNPS->pos==c->pos) allelecount[cellnumber][mySNPS->snpnumber][1]++;
          mySNPS=mySNPS->next;
        }
      }
      //Adds SNPS with alleles that are the ansestral type after the last SNP of derived type on the chromosome
      else if (c->type == 1) {
        while (mySNPS != NULL && mySNPS->pos < c->pos)
        {
          if (mySNPS->pos > start) allelecount[cellnumber][mySNPS->snpnumber][0]++;
          mySNPS=mySNPS->next;
        }
      }
      c = c->next;
    }
    G = G->nextChr;
  }
  char filename[30];
  FILE *SNPoutfile;
  printf("Printing alleles for cell %i\n",cellnumber);
  strcpy(filename, "true_snpAlleles_cancer_cell_");
  sprintf(filename, "%s%d",filename,cellnumber);
  if (NULL==(SNPoutfile=fopen(filename,"w")))
  {
    puts ("Cannot open true SNP allele count outfile!\n");
    exit(-1);
  }
  for (i=0; i<SNPNUMBER; i++) {
    //printf("SNP %i: %i %i\n",i,allelecount[cellnumber][i][0],allelecount[cellnumber][i][1]);
    fprintf(SNPoutfile,"%i\t%i\t%i\n",i,allelecount[cellnumber][i][0],allelecount[cellnumber][i][1]);
  }
  fclose(SNPoutfile);

  return temp;
}




void get_sim_params(char *paramname, struct param *parameters)
{
  FILE *paramfile;
  int i, even;
  double total;

  if (NULL==(paramfile=fopen(paramname,"r")))
  {
    puts ("Cannot open infile with parameters!");
    exit(-1);
  }
  fscanf(infile,"%i %i %lf %lf %lf %i\n",&parameters->numbins, &parameters->negbin_r, &parameters->number_reads, &parameters->readlength, &parameters->errorrate, &even);
  if (parameters->readlength> 0.1/(double)parameters->numbins)
  {
    printf("\n\nWith a read length of %lf and a bin length of %lf genomic units, the simulation assumptions are not met. ",parameters->readlength,1.0/(double)parameters->numbins);
    printf("The read length is assumed to be much smaller than the bin length.\n");\
      exit(-1);
  }
  if (parameters->numbins < 1 || parameters->number_reads <= 0.0) myerror(201);
  parameters->binrates = malloc(parameters->numbins*sizeof(double));
  if (even ==1) {
    for (i=0; i<parameters->numbins; i++)
      parameters->binrates[i]=1.0;
  }
  else {
    for (i=0; i<parameters->numbins; i++)
      fscanf(paramfile,"%lf",&parameters->binrates[i]);
    total=0.0;
    for (i=0; i<parameters->numbins; i++)
      total += parameters->binrates[i];
    for (i=0; i<parameters->numbins; i++)
      parameters->binrates[i]= (double)(parameters->numbins)*parameters->binrates[i]/total;
  }
  fclose(paramfile);
}

/*
   Notes on the probability that a read will cover a certain SNP:
   Let D be the total amount of DNA in the bin.
   let RL be the length of the read in the same units as D
   Let C be the CN state at the SNP
   Then the probability that a read from the bin hits the SNP is approx.
   C*RL/D
   This does not account for edge effects but is an infinitesimal approximation, i.e. assumes read length is infinitesimally small compared to bin length
   */


int globalcounter;
struct Genome* sim_and_print_binnumbers(int cellnumber, int cancer, struct Genome* G, struct param parameters, int ***allelecount, struct SNP* snps)
{

  double *bins, total, step, scale, m, p, AlleleFreq, af, tot, **SNPprob, **snppos;
  int i, j, k, x, D, A, s, snpCounter, *SNPcount, **snpnumber;
  char filename[30];
  FILE *outfile, *SNPoutfile;
  struct Genome* temp = G;
  struct SNP *snptemp, *snpstart;

  //all of this code in the beginning of the function is just to organize info about the SNPs in convenient data strcutures to facilitate the simulations

  step = GENOMELENGTH/parameters.numbins;
  bins=malloc(parameters.numbins*sizeof(double));
  for (i=0; i<parameters.numbins; i++)
    bins[i]=0;
  total = count_in_bin(G, parameters, bins);
  k=0;
  if (cancer==1){
    snpstart=snps;
    snps= snps->next;//first SNP is a dummy SNP at position 0
    SNPcount=malloc(parameters.numbins*sizeof(int));
    snpnumber=malloc(parameters.numbins*sizeof(int*));
    SNPprob=malloc(parameters.numbins*sizeof(double*));
    snppos=malloc(parameters.numbins*sizeof(double*));
    for (i=0; i<parameters.numbins; i++){
      SNPcount[i]=0;
      snptemp = snps;
      while (snps != NULL && snps->pos <(double)i*step){
        SNPcount[i]++;
        snps = snps->next;
      }
      snps = snptemp;
      SNPprob[i]=malloc(SNPcount[i]*sizeof(double));
      snpnumber[i]=malloc(SNPcount[i]*sizeof(int));
      snppos[i]=malloc(SNPcount[i]*sizeof(double*));
      j=0;
      // printf("AE %i %lf\n",i,snps->pos);
      while (snps != NULL && snps->pos <(double)i*step){
        // printf("Found snp %i in bin %i (%lf,%lf,%lf)\n",k,i,(double)(allelecount[cellnumber][k][0]+allelecount[cellnumber][k][1]),parameters.readlength,bins[i]);
        SNPprob[i][j] = (double)(allelecount[cellnumber][k][0]+allelecount[cellnumber][k][1])*parameters.readlength/bins[i];
        snpnumber[i][j]=k;
        snppos[i][j]=snps->pos ;
        j++;
        k++;
        snps = snps->next;
      }
    }
  }
  //prints out true mean ploidy in bin
  if (cancer==1) strcpy(filename, "true_copyNumber_cancer_cell_");
  else strcpy(filename,"true_copyNumber_healthy_cell_");
  sprintf(filename, "%s%d",filename,cellnumber);
  if (NULL==(outfile=fopen(filename,"w")))
  {
    puts ("Cannot open  read count outfile!\n");
    exit(-1);
  }
  if (cancer==1) {
    strcpy(filename, "simu_snpAlleles_cancer_cell_");
    sprintf(filename, "%s%d",filename,cellnumber);
    if (NULL==(SNPoutfile=fopen(filename,"w")))
    {
      puts ("Cannot open simu SNP allele read count outfile!\n");
      exit(-1);
    }
    snpCounter=0;
  }
  for (i=0; i<parameters.numbins; i++)
    fprintf(outfile,"chr1\t%.3lf\t%.3lf\t%.3lf\n",(double)i*step, (double)(i+1)*step, bins[i]/step);
  fclose(outfile);
  //simulates and prints out simulated cells
  if (cancer==1) strcpy(filename,"simu_readDepth_cancer_cell_");
  else strcpy(filename,"simu_readDepth_healthy_cell_");
  sprintf(filename, "%s%d",filename, cellnumber);
  if (NULL==(outfile=fopen(filename,"w")))
  {
    puts ("Cannot open outfile!\n");
    exit(-1);
  }
  scale = parameters.number_reads/total; //NOTICE I THINK SCALE CAN BE REMOVED
  tot=0;  //this step standardizes the expected number of reads to accomodate variation among bins
  for (i=0; i<parameters.numbins; i++)
  {
    m = parameters.binrates[i]*bins[i]*scale;
    tot+=m;
  }
  scale = scale*parameters.number_reads/tot;
  tot=0.0;  //below is where the simulations are actually happening
  for (i=0; i<parameters.numbins; i++)
  {
    m = parameters.binrates[i]*bins[i]*scale;
    p =  1.0 - m/((double)parameters.negbin_r + m);//notice p is 1-p from definition on wikipedia
    if (m>0.0){
      k=ignnbn(parameters.negbin_r, p);
      if (cancer==1){
        for (j=0; j<SNPcount[i]; j++)
        {
          D=A=0;
          x=binomialdev(SNPprob[i][j],k);
          //printf("%i %i %lf\n",x,k,SNPprob[i][j]);
          if (x>0) {
            af = (double)allelecount[cellnumber][snpnumber[i][j]][0]*(1.0-parameters.errorrate)+(double)allelecount[cellnumber][snpnumber[i][j]][1]*parameters.errorrate;
            AlleleFreq = af/((double)allelecount[cellnumber][snpnumber[i][j]][1]+af);
            //printf("AlleleFreq: %lf\n",AlleleFreq);
            A = binomialdev(AlleleFreq,x);
            D = x - A;
            //printf("SNP %i in bin %i: %i+%i=%i (snpprob=%lf, readdepth=%i)\n",j,i,A,D,x,SNPprob[i][j],k);
          }
          fprintf(SNPoutfile,"%i\t%i\t%lf\t%i\t%i\n",snpCounter,i,snppos[i][j],A,D);
          snpCounter++;
        }
      }
      fprintf(outfile,"chr1\t%.3lf\t%.3lf\t%i\n",(double)i*step, (double)(i+1)*step,k);
    }
    else if (cancer==1) {
      for (j=0; j<SNPcount[i]; j++) {
        fprintf(SNPoutfile,"%i\t%i\t%lf\t0\t0\n",snpCounter,i,snppos[i][j]);
        snpCounter++;
      }
      fprintf(outfile,"chr1\t%.3lf\t%.3lf\t0\n",(double)i*step, (double)(i+1)*step);
    }
  }
  if(cancer==1) {
    fclose(SNPoutfile);
  }
  fclose(outfile);

  if (cancer==1){
    for (i=0; i<parameters.numbins; i++){
      free(SNPprob[i]);
      free(snppos[i]);
      free(snpnumber[i]);
    }
    free(SNPprob);
    free(snppos);
    free(snpnumber);
    free(SNPcount);
    snps=snpstart;
  }
  free(bins);
  return(temp);
}

void sim_andprint_data(struct param parameters, int ***allelecount, struct SNP* snps, int cancer)
{
  int i;
  struct SNP* snptemp;
  snptemp = snps;
  for (i=0; i<numleafs; i++){
    sim_and_print_binnumbers(i, cancer, AllGenomes[i], parameters, allelecount, snps);
    snps=snptemp;
  }

}

//------------------------Tree based copy-number simulation algorithms-----------------------------



/*this algorithm recurses through the tree from the root to eventually generate the data in the leaf nodes*/
void recursesim(struct Genome* G, struct SNP* snps, int node, double lambda_de, double lambda_am, double length_de, double length_am, double theta, struct param parameters)
{
  //printf("Simulating branch of length %lf\n",tree[node].bl);
  G=simulate_edge(G, lambda_de,  lambda_am,  length_de,  length_am,  theta, tree[node].bl);
  if (node != -1 && tree[node].child[0]>-1)
  {
    struct Genome* G2 = (struct Genome*)malloc(sizeof(struct Genome));
    G2 = copyG(G);
    //test_genome(G2);
    recursesim(G, snps, tree[node].child[0], lambda_de,  lambda_am,  length_de,  length_am, theta, parameters);
    recursesim(G2, snps, tree[node].child[1], lambda_de,  lambda_am,  length_de,  length_am, theta, parameters);
  }
  else
  {
    /*printf("Genome in node %i:\n",node);
      printGenome(G);*/
    AllGenomes[globalcounter]=G;
    //    G=sim_and_print_binnumbers(globalcounter, 1, G, parameters);
    snps = find_SNPS(globalcounter, G, snps);
      /*printf("Current SNP set:\n");
        screen_print_SNPs(snps);
        printf("\n");*/
    /*print_SNPS(globalcounter, G);*/
    printf("Saving genome in node %i, cellNum %i\n",node, globalcounter); // Wed 26 Jan 2022 01:20:27 PM PST added
    globalcounter++;
    //delete_genome(G);
  }
}

/*core simulation algorithm*/
void sim_single_cell(double lambda_de, double lambda_am, double length_de, double length_am, double theta, double T, struct param parameters, int cancer)

{
  struct Genome* G;
  struct SNP* snps;
  int i, j;
  int ***allelecount;

  G = (struct Genome*)malloc(sizeof(struct Genome));
  makeDiploidGenomeNoGaps(G, GENOMELENGTH);

  if (cancer==1){
    snps  = (struct SNP*)malloc(sizeof(struct SNP));
    snps->pos = 0.0;
    snps->next = NULL;
    if (numleafs == 1) G=simulate_edge(G, lambda_de,  lambda_am,  length_de,  length_am,  theta, T);
    else recursesim(G, snps, root, lambda_de,  lambda_am,  length_de,  length_am, theta, parameters);
    allelecount = malloc(numleafs*(sizeof(int **)));
    for (i=0; i<numleafs; i++){
      allelecount[i] = malloc(SNPNUMBER*(sizeof(int *)));
      for (j=0; j<SNPNUMBER; j++){
        allelecount[i][j]= malloc(2*(sizeof(int)));
        allelecount[i][j][0]=allelecount[i][j][1]=0;
      }
    }
    snps=enumrate_SNPS(snps);
    /* printf("Current SNP set:\n");
       screen_print_SNPs(snps);*/
    for (i=0; i<numleafs; i++)
      snps=get_allelel_counts_in_bins(i, AllGenomes[i], snps, parameters.numbins, allelecount);
    sim_andprint_data(parameters, allelecount, snps, cancer);
    for (i=0; i<numleafs; i++){
      for (j=0; j<SNPNUMBER; j++)
        free(allelecount[i][j]);
      free(allelecount[i]);
    }
    free(allelecount);
    free(snps);
  }
  else {
    myerror(12);
    exit(-1);
  }
  free(G);

}


//-------------------FUNCTIONS FOR READING OR SIMULATING TREES-------------------------


void printtree()
{
  int i;

  for (i=0; i<2*numleafs-1; i++)
    printf("Node %i: children: (%i %i), parent: %i, bl: %.4lf\n",i,tree[i].child[0],tree[i].child[1],tree[i].parent,tree[i].bl);
}

void inittree()
{
  int i;

  for (i=0; i<numleafs; i++)
    tree[i].child[0]=tree[i].child[1]=-1;
  for (i=0; i<numleafs*2-1; i++)
    tree[i].bl=0.0;
}

// picks a time from the distribution of times to the next event with given rate
//algorithm from Slatkin and Hudson for exponential population growth
double slathud(double alpha, double rate, double thau)
{

  if (rate==0.0) myerror(200);
  if (alpha == 0.0) return expo(1.0);
  else return log(1.0 - alpha*exp(-thau)*(log(uniform())/rate))/alpha;
}

// Hudson's coalescence simulation algorithm
// returns the tmrca
double sim_coal_tree(double rootlength, double alpha)
{
  int i, j, pick, *list;
  double t, curt=0;

  printf("Simulating coalescence tree with %i leaf nodes\n",numleafs);
  inittree();
  list = malloc(numleafs*sizeof(int));
  for (i=0; i<numleafs; i++)
    list[i]=i;
  for (i=numleafs; i>1; i--){
    t = slathud(alpha, (double)(i*(i-1))/2, curt);
    curt += t;
    for (j=0; j<i; j++)
      tree[list[j]].bl += t;
    pick=i*uniform();
    tree[list[pick]].parent=2*numleafs-i;
    tree[2*numleafs-i].child[0]= list[pick];
    list[pick]=list[i-1];
    pick=(i-1)*uniform();
    tree[list[pick]].parent = 2*numleafs-i;
    tree[2*numleafs-i].child[1]= list[pick];
    list[pick]=2*numleafs-i;
  }
  root = 2*numleafs-2;
  tree[root].bl = rootlength;
  tree[root].parent = -1;
  free(list);
  //  printf("tmrca: %lf\n",t);
  return t;
}


/*allocates memory for the tree and reads it from the infile. Special code for the case of n=2*/
void setuptree(double T, int fromfile)
{
  double alpha;

  tree=malloc((numleafs*2-1)*(sizeof(struct treenode)));
  if (fromfile==1){
    if (numleafs == 2) //this is a bit awkward but the algorithm for reading the tree can't deal with >3 leaf nodes
    {
      tree[0].child[0]=1; tree[0].child[1]=2; tree[1].child[0]=tree[1].child[1]=tree[2].child[0]=tree[2].child[1]=-1;
      fscanf(infile,"%lf %lf\n",&tree[1].bl,&tree[2].bl);
      root=0;
    }
    else root=getclade()-1;
    if (root<0)
    {
      printf("Error reading tree!");
      exit(-1);
    }
    else tree[root].bl = T;
    printf("Read tree\n");
  }
  else {
    fscanf(infile,"%lf",&alpha);
    sim_coal_tree(T, alpha);
  }
}


//-------------------BIN BASED MODEL------------------------------------




void print_binnumbers(int cellnumber, int cancer, int *G, struct param parameters, double L)
{
  double tot, step, scale, m, p;
  char filename[30];
  int i;
  FILE *outfile;



  //prints out true mean ploidy in bin
  if (cancer==1) strcpy(filename, "true_copyNumber_cancer_cell_");
  else strcpy(filename,"true_copyNumber_healthy_cell_");
  sprintf(filename, "%s%d",filename,cellnumber);
  if (NULL==(outfile=fopen(filename,"w")))
  {
    puts ("Cannot open outfile!\n");
    exit(-1);
  }
  step = L/parameters.numbins;
  for (i=0; i<parameters.numbins; i++)
    fprintf(outfile,"chr1\t%.3lf\t%.3lf\t%.3lf\n",(double)i*step, (double)(i+1)*step, (double)G[i]);
  fclose(outfile);

  //simulates and prints out simulated cells
  if (cancer==1) strcpy(filename,"simu_readDepth_cancer_cell_");
  else strcpy(filename,"simu_readDepth_healthy_cell_");
  sprintf(filename, "%s%d",filename, cellnumber);
  if (NULL==(outfile=fopen(filename,"w")))
  {
    puts ("Cannot open outfile!\n");
    exit(-1);
  }

  tot=0;  //this step standardizes the expected number of reads to accomodate variation among bins
  for (i=0; i<parameters.numbins; i++)
    tot += parameters.binrates[i]*G[i];
  scale = parameters.number_reads/tot;
  tot=0.0;  //this is very the simulations are actually happening
  for (i=0; i<parameters.numbins; i++)
  {
    m = parameters.binrates[i]*G[i]*scale;
    p =  1.0 - m/(parameters.negbin_r + m);//notice p is 1-p from definition on wikipedia
    if (m>0.0)
      fprintf(outfile,"chr1\t%.3lf\t%.3lf\t%i\n",(double)i*step, (double)(i+1)*step,ignnbn(parameters.negbin_r, p));
    else
      fprintf(outfile,"chr1\t%.3lf\t%.3lf\t0\n",(double)i*step, (double)(i+1)*step);
  }
  fclose(outfile);
}

void mutation(int start, int end, int size, int max, int *G)
{
  int i, d;

  printf("Mutation from %i to %i of size %i\n",start,end,size);

  for (i=start; i<end+1; i++)
  {
    d=G[i]+size;
    if (d <= 0) G[i]=0;
    else if (d >= max) G[i] =  max;
    else G[i] = d;
  }
}


void mutate(double **cumurates, int *G, double p, int max, int numbins, double rateinternal, double rate)
{
  int i, j, end, start;
  double r;


  if (uniform()< rateinternal/rate) start = (int)(uniform()*(numbins-1))+1;
  else start = -1;
  end = start;
  while (end < numbins-1 && uniform()>p) end++;
  if (start==-1  && end -1);
  else {
    if (start==-1) start = 0;
    r = uniform();
    i=0;
    j= G[start];
    if (cumurates[j][max]==0.0) i=j;
    else {
      while (r>cumurates[j][i]) i++;
      printf("r: %lf, cumurates: %lf, j: %i, i: %i\n",r,cumurates[j][i],j,i);
      mutation(start, end, i-j, max, G);
    }
  }
}

void simulate_BINedge(double **ratematrix, double **cumurates, int *G, double p, int max, int numbins, double T)
{
  int i, j, k;
  double *binrates, rateinternal, rate, t = 0.0;

  binrates = malloc(numbins*(sizeof(double)));
  rateinternal = 0;
  for (i=0; i<numbins; i++)
  {
    binrates[i]=0;
    k = G[i];
    for (j=0; j<=max; j++)
    {
      if (k != j) binrates[i] += ratematrix[k][j];
    }
    rateinternal += binrates[i];
  }
  rate = rateinternal + rateinternal/(p*(double)numbins);//this is the trick to ensure a uniform rate across the genome
  t=0;
  do {
    t += expo(rate);
    printf("Time %lf, rates: %lf %lf, T: %lf\n",t,rate,rateinternal,T);
    if (t<T) {
      mutate(cumurates, G,  p,  max,  numbins,  rateinternal,  rate);
      //instead of this redundant it would be better with some code tghat faster calculated the change in the rate
      rateinternal = 0;
      for (i=0; i<numbins; i++)
      {
        binrates[i]=0;
        k = G[i];
        for (j=0; j<=max; j++)
        {
          if (k != j) binrates[i] += ratematrix[k][j];
        }
        rateinternal += binrates[i];
      }
      rate = rateinternal + rateinternal/(p*(double)numbins);//this is the trick to ensure a uniform rate across the genome
    }
  } while (t<T);
  free(binrates); //maybe make this a global or pass it around with pointers
}

void setup_healthy(int numbins, int G[])
{
  int i;

  for (i=0; i<numbins; i++)
    G[i]=2;
}



void copyBIN(int *G, int *G2, int numbins)
{
  int i;

  for (i=0; i<numbins; i++)
    G2[i] = G[i];
}

//this algorithm recurses through the tree from the root to eventually generate the data in the leaf nodes for binned data
void recursesimBIN(int *G, int node, double **cumumatrix, double **ratematrix, double p, struct param parameters, int max, double L)
{

  simulate_BINedge(ratematrix, cumumatrix, G, p, max, parameters.numbins, tree[node].bl);
  if (tree[node].child[0]>-1)
  {
    int *G2 = malloc(parameters.numbins*sizeof(int));
    copyBIN(G, G2, parameters.numbins);
    recursesimBIN(G, tree[node].child[0], cumumatrix, ratematrix, p, parameters, max, L);
    recursesimBIN(G2, tree[node].child[1], cumumatrix, ratematrix, p, parameters, max, L);
  }
  else
  {
    printf("Saving node %i, cellNum %i:\n",node, globalcounter);
    print_binnumbers(globalcounter, 1, G, parameters, L);
    globalcounter++;
    free(G);
  }
}

//core simulation algorithm
void sim_single_cell_bin(double T, double **ratematrix, double p, struct param parameters, int cancer, int max, double L)

{
  int i, j, *G;
  double **cumurates;

  cumurates = malloc((max+1)*(sizeof(double *)));
  for (i=0; i<=max; i++)
    cumurates[i]= malloc((max+1)*(sizeof(double)));
  for (i=0; i<=max; i++)
  {
    cumurates[i][0]=ratematrix[i][0];
    ratematrix[i][0]=ratematrix[i][0]/(double)parameters.numbins;
    for (j=1; j<=max; j++){
      cumurates[i][j]= cumurates[i][j-1];
      if (i!=j) cumurates[i][j] += ratematrix[i][j];
      ratematrix[i][j]=ratematrix[i][j]/(double)parameters.numbins;
    }
    for (j=0; j<=max; j++){
      if (cumurates[i][max]>0.0)
        cumurates[i][j]=cumurates[i][j]/cumurates[i][max];
    }

  }

  G = malloc(parameters.numbins*sizeof(int));
  setup_healthy(parameters.numbins,  G);
  if (numleafs == 1)
  {
    simulate_BINedge(ratematrix, cumurates, G, p,  max, parameters.numbins, T);
    print_binnumbers(globalcounter, cancer, G, parameters, L);
    globalcounter++;
    free(G);
  }
  else recursesimBIN(G, root, cumurates, ratematrix, p, parameters, max, L);
  for (i=0; i<=max; i++)
    free(cumurates[i]);
  free(cumurates);
}

void checkratematrix(double **ratematrix, int max)
{
  int i, j;
  double s, e=0.000000001;

  for (i=0; i<=max; i++)
  {
    s=0.0;
    for (j=0; j<=max; j++)
    {
      s+=ratematrix[i][j];
    }
    if (s<0-e || s > e)
    {
      printf("Rate matrix row sums are non-zero! Can't have that.");
      exit(-1);
    }
  }
}


//------------------------------------------------------------------------------

/*Reads two files

  First file includes paramaters for cell evolution simulation and possibly a tree
  infile format:
  [if bin based model (1: bin based)] [Genomelength] [if tree should be simulated or from file (1: from file)] [number of healthy cells to siumulate]
  the next lines include parameters for the model of evolution If the bin based model is NOT chosen the line has the format
  [rate of deletion] [rate of amplifications] [mean deletion length] [mean amplification length] [theta=mutation rate per time unit per chromosome length] [per SNP allele error rate] [number of cells in tree]

  If the bin based model has been chosen:
  [max ploidy] [parameter p] (p is the probability of NOT extending to the next bin)
  [first row of rate matrix]
  [second row of rate matrix]
  ...
  [last row of rate matrix] (this is a matrix of dimension (maxploidy+1)*(maxploidy+1))


  The follows specification of the tree and branchlengths:
  [length of edge leading to root of tree = length of simulation time for single cell]
  [newick format treefile with branchlengths. if there are exactly two cells instead of a tree just give the branclengths of the two bracnhes]

  Second file includes parameters for simulating read counts in bins.  See paramfile.txt and paramfile2.txt for eaxmples


  Output files:
  true_copyNumber_cancer_cell_i: true copy number in each bin [chr start end CN]
  simu_readDepth_cancer_cell_i: simulated read depth in each bin [chr start end readDepth]

  true_snpAlleles_cancer_cell_i: true allele freq for each snp [idx #ancestralAlleles #derivedAlleles]
  simu_snpAlleles_cancer_cell_i: simulated read depth for each snp [idx bin# posInSimGenome #ancReads #derReads]

  true_copyNumber_healthy_cell_i: true copy number in each bin (should be all 2 for healthy diploid cells) [chr start end CN]
  simu_readDepth_healthy_cell_i: simulated read depth in each bin [chr start end readDepth]
  */


//compile using gcc CNV.c ranlib.c rnglib.c -o CNV -lm
int main(int argc, char *argv[])
{

  double lambda_deletion, lambda_amplification, length_deletion, length_amplification, theta, p, T, **ratematrix;
  int i, j, fromfile, numberhealthy, binbased, max;
  char inname[30], paramname[30];
  struct param simparam;
  globalcounter=0;

  //SetSeed(23);
  if (argc<3) {printf("Specify name of two infiles: infile paramfile [seed]\n"); exit(-1);}
  sprintf(inname, "%s", argv[1]);
  sprintf(paramname, "%s", argv[2]);
  if(argc == 4) {
    SetSeed(atoi(argv[3]));
  } else {
    SetSeed(time(NULL));
  }
  if (NULL==(infile=fopen(inname,"r")))
  {
    puts ("Cannot open infile with!");
    exit(-1);
  }
  fscanf(infile,"%i %lf %i %i %i\n",&binbased, &GENOMELENGTH, &numleafs, &fromfile, &numberhealthy);
  if (binbased)
  {
    fscanf(infile,"%i %lf\n", &max, &p);
    ratematrix = malloc((max+1)*(sizeof(double *)));
    for (i=0; i<=max; i++)
      ratematrix[i]= malloc((max+1)*(sizeof(double)));
    for (i=0; i<=max; i++)
      for  (j=0; j<=max; j++)
        fscanf(infile,"%lf\n",&ratematrix[i][j]);
    checkratematrix(ratematrix, max);
  }
  else
    fscanf(infile,"%lf %lf %lf %lf %lf\n",&lambda_deletion, &lambda_amplification, &length_deletion, &length_amplification, &theta);
  fscanf(infile,"%lf\n",&T);
  if (numleafs > 1) setuptree(T, fromfile);
  //printtree();
  fclose(infile);
  get_sim_params(paramname, &simparam);
  if (binbased)//cannot simulate SNPS
  {
    printf("Now simulating %i cancer cell(s) using bin based  model\n",numleafs);
    sim_single_cell_bin(T, ratematrix,  p, simparam, 1, max, GENOMELENGTH);
    if (numberhealthy >0)
    {
      globalcounter=0;
      numleafs=1;
      printf("Now simulating %i healthy cell(s)\n",numberhealthy);
      for (i=0; i<numberhealthy; i++)
        sim_single_cell_bin(0.0, ratematrix, p, simparam, 0, max, GENOMELENGTH);
    }
  }
  else {
    SNPNUMBER=0;
    printf("Now simulating %i cancer cell(s) using line segment model\n",numleafs);
    /*allocate array of genomes*/
    AllGenomes = malloc((numleafs*2-1)*sizeof(struct Genome));
    sim_single_cell(lambda_deletion, lambda_amplification, length_deletion, length_amplification, theta, T, simparam, 1);
    printtree();
    if (numberhealthy >0)
    {
      struct Genome* G;
      struct SNP* snps;
      int ***allelecount;
      globalcounter=0;
      numleafs=1;
      G = (struct Genome*)malloc(sizeof(struct Genome));
      makeDiploidGenomeNoGaps(G, GENOMELENGTH); // use this version for regular simulations
      //makeDiploidGenomeNoGapsDuplicated(G, GENOMELENGTH); // use this version for a whole genome duplication (ie starts off with 4 chr instead of 2)
      printf("Now simulating %i healthy cell(s)\n",numberhealthy);
      for (globalcounter=0; globalcounter<numberhealthy; globalcounter++){
        sim_and_print_binnumbers(globalcounter, 0, G, simparam, allelecount, snps);
      }
    }
  }
}


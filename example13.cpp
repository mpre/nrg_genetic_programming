#include <iostream>
#include <cmath>
#include <fstream>
#include <iostream>

#include <ga/ga.h>
#include <ga/GASelector.h>
#include <ga/GATree.C>
#include <ga/GATreeGenome.C>

#define MAX_DEPTH 5
#define MIN_DEPTH 1
#define FUNC_N 5
#define TERM_N 8
#define POP_SIZE 100

#define SOGLIA 20000

int initialized_genomes =0;
int depth_of_this_bucket=1;
class node_content
{
public:
  node_content()
  {
    terminal = true;
    type = 0;
  }
  node_content(bool ter, int t)
  {
    terminal = ter;
    type = t;
  }
  node_content(const node_content& n)
  {
    terminal = n.terminal;
    type = n.type;
  }
  node_content& operator=(const node_content& n) {
    terminal = n.terminal;
    type = n.type;
    return *this;
  }
  ~node_content() {}
  /* indica se il nodo è terminale o funzionale */
  bool terminal;
  int type;
};

template class GATreeGenome<node_content>;
template class GATree<node_content>;

float** values;
int nrows;
int nvars;

void init_values()
{
  std::ifstream f;
  f.open("energia.txt");
  f >> nrows;
  f >> nvars;
  values = new float*[nrows];
  for (int i = 0; i < nrows; i++) {
    values[i] = new float[nvars];
    for (int j = 0; j < nvars; j++) {
      f >> values[i][j];
    }
  }
  f.close();
}

float eval(GATreeGenome<node_content>& tree, int n)
{
  if (tree.current()->terminal) {
    return values[n][tree.current()->type];
  } else {
    int type = tree.current()->type;
    tree.child();
    double tmp1 = eval(tree,n);
    tree.next();
    double tmp2 = eval(tree,n);
    double tmp3 = 1;
    if (type == 4) {  // IF
      tree.next();
      tmp3 = eval(tree, n);
    }
    tree.parent();
    switch (type) {
    case 0:
      return tmp1 + tmp2; break;
    case 1:
      return tmp1 - tmp2; break;
    case 2:
      return tmp1 * tmp2; break;
    case 3:
      return (tmp2!=0?tmp1/tmp2:0); break;
    case 4:
      if (tmp1 > SOGLIA)
	return tmp2;
      else
	return tmp3;
   }
    return -1;
  }
}

void print_tree(const GAGenome& g);

float objective(GAGenome& g)
{
  GATreeGenome<node_content>& tree = (GATreeGenome<node_content> &)g;
  double score = 0;
  float modificator = 1;
  tree.root();
  for (int i = 0; i < nrows; i++)
    score += std::pow(values[i][nvars-1] - eval(tree,i),2);
  if (tree.depth() > MAX_DEPTH)
    modificator = 4 * (float)tree.depth() / (float)MAX_DEPTH;
  if (tree.depth() <= MIN_DEPTH)
    modificator = 3;
  return std::sqrt(score/nrows) * modificator ;
}

void init_random(GATreeGenome<node_content>& tree, int depth)
{
  if ((depth == 0) || GAFlipCoin(0.5)) {
    node_content n(true,GARandomInt(0,TERM_N));
    tree.insert(n,GATreeBASE::BELOW);
  } else {
    node_content n(false,GARandomInt(0,FUNC_N));
    tree.insert(n,GATreeBASE::BELOW);
    int arity = n.type==4?3:2;
    for (int i = 0; i < arity; i++)
      init_random(tree,depth-1);
  }
  tree.parent(); //riporto l'iteratore al posto corretto.
}

/* Inizializzazione grow */
void init_tree_grow(GAGenome& g)
{
  GATreeGenome<node_content>& tree = (GATreeGenome<node_content> &)g;
  tree.root();
  tree.destroy();
  node_content n(false,GARandomInt(0,FUNC_N));
  tree.insert(n,GATreeBASE::ROOT);
  int arity = n.type==4?3:2;
  for (int i = 0; i < arity; i++)
    init_random(tree,depth_of_this_bucket);
}

void init_full(GATreeGenome<node_content>& tree, int depth) {
  if(depth==0) {
    node_content n(true, GARandomInt(0, TERM_N));
    tree.insert(n, GATreeBASE::BELOW);
  } else {
    node_content n(false, GARandomInt(0, FUNC_N));
    tree.insert(n, GATreeBASE::BELOW);
    int arity = n.type==4?3:2;
    for (int i=0; i<arity; ++i)
      init_full(tree, depth-1);
  }
  tree.parent();
}

/* Inizializzazione full */
void init_tree_full(GAGenome& g)
{
  GATreeGenome<node_content>& tree = (GATreeGenome<node_content>&) g;
  tree.root();
  tree.destroy();
  node_content n(false, GARandomInt(0, FUNC_N));
  tree.insert(n, GATreeBASE::ROOT);
  int arity = tree.current()->type==4?3:2;
  for(int i=0; i<arity; ++i) {
    init_full(tree, depth_of_this_bucket);
  }
}

/* Inizializzazione ramped half-half */
void init_ramped_half_half(GAGenome& g) {
  
  int elements_in_a_bucket = POP_SIZE / (MAX_DEPTH*2);
  int full_buckets = initialized_genomes / elements_in_a_bucket;
  int elements_in_this_bucket = initialized_genomes - elements_in_a_bucket*full_buckets;
  depth_of_this_bucket = full_buckets + 1;
  std::cout << "Elementi in questo bucket " << elements_in_this_bucket << " su " << elements_in_a_bucket << std::endl;
  bool grow_or_full;
  
  std::cout << "Inizializzato un ";
  if (elements_in_this_bucket < elements_in_a_bucket/2)
    grow_or_full = true;
  else
    grow_or_full = false;
  if(grow_or_full) {
    init_tree_grow(g);
    std::cout << "grow " << std::endl;
  }
  else { 
    init_tree_full(g);
    std::cout << "full " << std::endl;
  }
  std::cout << "Depth : " << depth_of_this_bucket << std::endl;
  initialized_genomes +=1;
  std::cout << "Inizializzati : " << initialized_genomes << std::endl;
  print_tree(g);
  std::cout <<std::endl << std::endl;
}

void print_node_content(const GAGenome& g)
{
  GATreeGenome<node_content>& tree = (GATreeGenome<node_content> &)g;
  if (tree.current()->terminal) {
    std::cout << "x" << tree.current()->type;
  } else {
    std::cout << "(";
    switch (tree.current()->type) {
    case 0:
      std::cout << "+ "; break;
    case 1:
      std::cout << "- "; break;
    case 2:
      std::cout << "* "; break;
    case 3:
      std::cout << "/ "; break;
    case 4:
      std::cout << "if (> ";      
    }
  }
}

void print_tree(const GAGenome& g)
{
  GATreeGenome<node_content>& tree = (GATreeGenome<node_content> &)g;
  print_node_content(tree);
  if (!tree.current()->terminal) {
    int if_test = 0;
    if (tree.current()->type == 4)
      if_test = 1;
    tree.child();
    print_tree(tree);
    if(if_test) { 
      std::cout << " " << SOGLIA << ")";
    }
    std::cout << " ";
    tree.next();
    print_tree(tree);
    if(if_test){
      std::cout << " ";
      tree.next();
      print_tree(tree);
    }
    tree.parent();
    std::cout << ")";
  }
}

int main()
{
  init_values();
  GARandomSeed();
  GATreeGenome<node_content> genome(objective);
  genome.initializer(init_ramped_half_half);
  genome.crossover(GATreeGenome<node_content>::OnePointCrossover);
  genome.mutator(GATreeGenome<node_content>::SwapSubtreeMutator);
  GASimpleGA ga(genome);
  ga.initialize();
  ga.minimize();
  ga.populationSize(POP_SIZE);
  ga.pMutation(0.1);
  ga.pCrossover(0.9);
  ga.elitist(gaTrue);
  for (int i = 0; initialized_genomes * i < 10000; i++) {
    ++ga;
    print_tree(ga.statistics().bestIndividual());
    std::cout << std::endl;
    std::cout << objective(const_cast<GAGenome&>(ga.statistics().bestIndividual())) << std::endl;
    //std::cout << std::endl;
  }
}

#include <iostream>
#include <cmath>
#include <fstream>

#include <ga/ga.h>
#include <ga/GASelector.h>
#include <ga/GATree.C>
#include <ga/GATreeGenome.C>

#define MAX_DEPTH 8
#define MIN_DEPTH 2
#define FUNC_N 5

#define SOGLIA 20000

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
    if (n <= nvars -2)
      return values[n][tree.current()->type];
    else{ 
      switch (n + 2 - nvars) {
      case 0:
	return 1000;
	break;
      case 1:
	return 2000;
	break;
      }
    }
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
    node_content n(true,GARandomInt(0,nvars-2));
    tree.insert(n,GATreeBASE::BELOW);
  } else {
    node_content n(false,GARandomInt(0,FUNC_N));
    tree.insert(n,GATreeBASE::BELOW);
    for (int i = 0; i < 2; i++)
      init_random(tree,depth-1);
  }
  tree.parent(); //riporto l'iteratore al posto corretto.
}

/* Inizializzazione random */
void init_tree(GAGenome& g)
{
  GATreeGenome<node_content>& tree = (GATreeGenome<node_content> &)g;
  tree.root();
  tree.destroy();
  node_content n(false,GARandomInt(0,FUNC_N));
  tree.insert(n,GATreeBASE::ROOT);
  for (int i = 0; i < 2; i++)
    init_random(tree,MAX_DEPTH-1);
}

void print_node_content(const GAGenome& g)
{
  GATreeGenome<node_content>& tree = (GATreeGenome<node_content> &)g;
  if (tree.current()->terminal) {
    if(tree.current()->type <= nvars - 2)
      std::cout << "x" << tree.current()->type;
    else 
      switch(tree.current()->type - nvars + 2) { 
      case 0:
	std::cout << "1000";
	break;
      case 1:
	std::cout << "2000";
	break;
      }
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
  genome.initializer(init_tree);
  genome.crossover(GATreeGenome<node_content>::OnePointCrossover);
  genome.mutator(GATreeGenome<node_content>::SwapSubtreeMutator);
  GASimpleGA ga(genome);
  ga.initialize();
  ga.minimize();
  ga.populationSize(100);
  ga.pMutation(0.1);
  ga.pCrossover(0.9);
  ga.elitist(gaTrue);
  for (int i = 0; i < 100; i++) {
    ++ga;
    print_tree(ga.statistics().bestIndividual());
    std::cout << std::endl;
    std::cout << objective(const_cast<GAGenome&>(ga.statistics().bestIndividual())) << std::endl;
    //std::cout << std::endl;
  }
}

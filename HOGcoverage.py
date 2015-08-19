# import external libraries
import random, re, dendropy, familyanalyzer as fa
import lxml.etree as etree
import sys, getopt

def main(argv):

    try:
        opts, args = getopt.getopt(argv,"f:")
    except getopt.GetoptError:
        print('Usage: countHOGs.py -f [PATH/TO/FILE] ')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-f"):
            XML = arg

    # parse OrthoXML file, compute taxonomy and apply to data
    op = fa.OrthoXMLParser(XML)

    tax = fa.TaxonomyFactory.newTaxonomy(op)
    op.augmentTaxonomyInfo(tax)
    root = tax.root
    # compute all events on all branches
    tax.get_histories(op)
    all_comparisons = tax.get_comparisons(op)
    gtt = fa.GeneTreeTracer(op, tax)
    toplevelog = op.getToplevelGroups()
    hoglist = [fa.GeneFamily(fam) for fam in toplevelog]
    # compute number of losses and number of leaves for each HOG
    hogcounts = dict()
    for hog in hoglist:
        hogcounts[hog.getFamId()] = dict()
        levels = hog.getLevels()

        # get HOG tree
        gene_tree = gtt.trace_gene_family(tax[levels[0]], hog, None, False)
        newick_long = gene_tree.write()

        #remove meta data
        newick_short = re.sub(r'\[.*?\]', '',newick_long.replace('"',"'"))+';'

        # count losses and leaves
        losses = newick_long.count('LOSS');
        tree = dendropy.Tree.get_from_string(newick_short, schema="newick")
        retains = len(list(tree.leaf_node_iter()))
        hogcounts[hog.getFamId()] = {"losses":losses, "leaves": retains}

    # Output
    print("HOG ID,#sequences,HOG coverage")
    for hog in hoglist:
        famId = hog.getFamId()
        d = 1.0*(hogcounts[famId]['leaves']-hogcounts[famId]['losses'])/hogcounts[famId]['leaves']
        print(famId + "," + str(len(hog.getMemberGenes())) + "," + str(d))


if __name__ == "__main__":
   main(sys.argv[1:])
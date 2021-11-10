#!/usr/bin/env python3
"""
compare_trees.py is a tool to compare a tree to a reference
tree and show the discrepancies between the two. If taxon
metadata is provided it shows the violation of monophyletic
topology too.

Author: Lenard Szantho <lenard@drenal.eu>
Version: 
    - 0.1 (2020-04-07)
    - 0.2 (2021-11-10) more robust leaf name matching, removed blue circles
"""

import os
import argparse
import csv
import itertools
import re

from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces
from colorhash import ColorHash


def main():
    """ Entrypoint of program """

    # set and parse arguments
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "tree", help="path to tree of question (newick or NGX newick format)",
    )
    parser.add_argument(
        "reftree",
        help="path to the reference tree to compare with (newick or NGX newick format)",
    )
    parser.add_argument(
        "outgroup",
        metavar="outgroup",
        nargs="+",
        help="List of species (leaves) to use as outgroup when rooting the trees",
    )
    parser.add_argument(
        "--taxa",
        dest="taxa",
        help='metadata file with species names and taxa (csv file with syntax: "species","taxa")',
    )
    parser.add_argument(
        "--taxatree",
        dest="taxatree",
        help="newick tree representing the supposed topology of the taxa",
    )
    parser.add_argument(
        "--outputprefix",
        dest="outputprefix",
        help="attach prefix to the output files (by default: tree.png, reftree.png, report.txt)",
        default="",
    )
    args = parser.parse_args()

    # read in trees
    tree = Tree(args.tree, format=1)
    reftree = Tree(args.reftree, format=1)

    # outgroup them
    if len(args.outgroup) > 1:
        tree_outgroup_ancestor = tree.get_common_ancestor(args.outgroup)
        reftree_outgroup_ancestor = reftree.get_common_ancestor(args.outgroup)
        tree.set_outgroup(tree_outgroup_ancestor)
        reftree.set_outgroup(reftree_outgroup_ancestor)
    else:
        tree.set_outgroup(args.outgroup[0])
        reftree.set_outgroup(args.outgroup[0])

    # calculate Robinson-Foulds distance between the two trees
    common_leaves = calculate_robinson_foulds(tree, reftree)

    # some basic comparison
    print(
        "Number of leaves in tree {}, reftree {}. Common leaves: {}".format(
            len(tree), len(reftree), len(common_leaves)
        )
    )

    # compare leaf pairs and mark discrepancies
    compare_leaf_pairs(tree, reftree)

    # -----------------------------
    # in case taxa file is supplied
    # -----------------------------
    species = set()
    taxa2species = dict()
    species2taxa = dict()
    logical_taxa = set()
    if args.taxa:
        # parse taxon data
        with open(args.taxa, mode="r") as taxacsv:
            entries = csv.reader(taxacsv, delimiter=",", quotechar='"')
            for row in entries:
                species.add(row[0])
                taxa2species.setdefault(row[1], []).append(row[0])
                species2taxa[row[0]] = row[1]

        # if species is a key in taxa, then it's not a species but a higher taxon
        for taxon in taxa2species:
            if taxon in species:
                logical_taxa.add(taxon)
                species.remove(taxon)

        # expand species names to their full length as given in the taxa csv
        restore_leaf_names(tree, species)
        restore_leaf_names(reftree, species)

        # the taxa dictionary may contain more species than the trees
        ignored_species_tree = search_nonpresent_species(tree, species)
        ignored_species_reftree = search_nonpresent_species(reftree, species)

        if ignored_species_tree:
            print(
                "Following species weren't found in the tree, they will be ignored: \n{}".format(
                    ", ".join(ignored_species_tree)
                )
            )

        if ignored_species_reftree:
            print(
                "Following species weren't found in the reference tree, they will be ignored: \n{}".format(
                    ", ".join(ignored_species_reftree)
                )
            )

        # add feature to leaves describing their taxon hierarchy
        identify_taxon_hierarchy(tree, species2taxa)
        identify_taxon_hierarchy(reftree, species2taxa)

        # check monophyly of taxa
        print("Checking monophyly of tree:\n")
        check_monophyly(tree, species, logical_taxa, taxa2species, ignored_species_tree)
        print("Checking monophyly of reference tree:\n")
        check_monophyly(
            reftree, species, logical_taxa, taxa2species, ignored_species_reftree
        )

    # -----------------------
    # if taxatree is supplied
    # -----------------------
    taxatree = Tree()
    if args.taxatree:
        taxatree = Tree(args.taxatree, format=1)

    # ----------
    # outputting
    # ----------
    output_tree(args.tree, tree, False, args.outputprefix)
    output_tree(args.reftree, reftree, True, args.outputprefix)

    # END


def output_tree(treefile, tree, ref=False, prefix=""):
    """
Output tree as a PNG file with proper formatting and
all the infromation gathered presented in a visual way
    """
    treefile_base = os.path.basename(treefile)

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = True
    ts.show_branch_support = False
    ts.show_scale = True
    ts.show_border = False
    ts.mode = "r"
    ts.min_leaf_separation = 1
    # this makes the tree too wide
    # ts.optimal_scale_level = "full"
    # ts.scale = 0.01
    #if ref:
    #    ts.title.add_face(
    #        faces.TextFace("Reference tree\n(file: {})".format(treefile_base)), column=0
    #    )
    #else:
    #    ts.title.add_face(
    #        faces.TextFace("Tree in question\n(file: {})".format(treefile_base)),
    #        column=0,
    #    )

    for node in tree.traverse("postorder"):
        if node.is_leaf():
            if "taxa" in node.features:
                c = ColorHash(node.taxa[0])
                node.add_feature("colour", c.hex)

    ts.layout_fn = ete3_layout
    if ref:
        tree.render(prefix + "reftree.png", w=1280, tree_style=ts)
    else:
        tree.render(prefix + "tree.png", w=1280, tree_style=ts)


def ete3_layout(node):
    nstyle = NodeStyle()

    # setting dots on every node
    nstyle["size"] = 0
    nstyle["shape"] = "circle"

    # discrepancy marker
    try:
        if node.discrepancy:
            nstyle["fgcolor"] = "red"
            nstyle["size"] = 10
        else:
            nstyle["fgcolor"] = "blue"
    except:
        pass

    if node.is_leaf():
        face_name = faces.AttrFace("name", fsize=10, fgcolor="#000000")
        faces.add_face_to_node(face_name, node, column=0)

        if "colour" in node.features:
            nstyle["bgcolor"] = node.colour

        if "taxa" in node.features:
            taxon = node.taxa[0]
            face_taxon = faces.TextFace(", " + taxon, fsize=10)
            faces.add_face_to_node(face_taxon, node, column=1)

        # face_clade_type = faces.AttrFace("clade_type", fsize=10, fgcolor="#000000")
        # faces.add_face_to_node(face_clade_type, node, column=3)
    else:
        if node.name != "NoName":
            face_name = faces.TextFace(node.name, fsize=8, fgcolor="#000000")
            faces.add_face_to_node(face_name, node, column=0, position="branch-bottom")

    node.set_style(nstyle)


def calculate_robinson_foulds(tree, reftree):
    """
Calculate Robinson-Foulds distance between the two trees

Parameters
----------
tree
    ete3 Tree instance
reftree
    ete3 Tree instance

Returns
-------
int : number of common leaves
    """
    (
        rf,
        max_rf,
        common_leaves,
        collapsed_leaves1,
        collapsed_leaves2,
        parts_t1,
        parts_t2,
    ) = reftree.robinson_foulds(tree)

    # print("sth1 {}".format(collapsed_leaves1))
    # print("sth2 {}".format(collapsed_leaves2))

    print("Robinson-Foulds distance is {} over a total of {}".format(rf, max_rf))
    print(
        "Partitions of tree not found in reference tree: {}".format(
            len(parts_t2) - len(parts_t1)
        )
    )
    print(
        "Partitions of reference tree not found in tree: {}".format(
            len(parts_t1) - len(parts_t2)
        )
    )

    return common_leaves


def compare_leaf_pairs(tree, reftree):
    """
Calculate Robinson-Foulds distance between the two trees

Parameters
----------
tree
    ete3 Tree instance
reftree
    ete3 Tree instance
    """
    # create all possible leaf pairs
    # leaves = reftree.get_leaf_names()

    # leaf_pairs = list(itertools.combinations(leaves, 2))

    tag_nodes(tree)
    tag_nodes(reftree)

    for node in tree.iter_descendants("postorder"):
        if not node.is_leaf():
            reftree_node = reftree.search_nodes(descendants=node.descendants)
            # reftree_node = search_polytomic_node(reftree, node.descendants)
            if not reftree_node:

                node.add_feature("discrepancy", True)
                iterate_subnodes(node, reftree)


def search_polytomic_node(reftree, descendants):
    """
Finds nodes who have the specified descentants, but allows
having more descendants than those listed.
    """
    for n in reftree.traverse("postorder"):
        all_found = True
        for l in descendants:
            if l not in n.descendants:
                all_found = False
        if all_found:
            return n
    return None


def iterate_subnodes(node, reftree):
    """
Iterate subnodes recursively until all discrepancies are found
    """
    for sub_node in node.iter_descendants("postorder"):
        if not sub_node.discrepancy:
            reftree_node = reftree.search_nodes(descendants=sub_node.descendants)
            # reftree_node = search_polytomic_node(reftree, node)

            if not reftree_node:
                sub_node.add_feature("discrepancy", True)
                # stop if it's a leaf, iterate further if internal
                if not sub_node.is_leaf():
                    iterate_subnodes(sub_node, reftree)


def tag_nodes(tree):
    """
Tag internal nodes by their childrens names

Parameters
----------
tree
    ete3 Tree instance
    """
    for node in tree.traverse("postorder"):
        # the discrepancy marker will tell if tree and referencetree doesn't match
        node.add_feature("discrepancy", False)

        # the descendants feature holds all descendants' name in a set
        if not node.is_leaf():
            descendants = set()
            for c in node.get_children():
                descendants |= c.descendants

            node.add_feature("descendants", descendants)
        #            print("New descendants: {}".format(descendants))
        else:
            node.add_feature("descendants", {node.name})


def restore_leaf_names(tree, species):
    """
The names of leaves may be cutted because
of the phylogenetic tools used,
this function restores the original names.

Parameters
----------
tree
    ete3 Tree instance
species : set
    list of species

    """
    leaf_name_length = max([len(leaf.name) for leaf in tree.get_leaves()])

    for leaf in tree.get_leaves():
        for spec in [s for s in species if len(s) > leaf_name_length]:
            if leaf.name == spec[:leaf_name_length]:
                # change the name to the full name
                leaf.name = spec


def search_nonpresent_species(tree, species):
    """
Collect the nonpresent species in the tree to a set

Parameters
----------
tree
    ete3 Tree instance
species : set
    list of species

    """
    to_ignore = set()
    leaves = tree.get_leaf_names()
    for spec in species:
        if spec not in leaves:
            to_ignore.add(spec)

    return to_ignore


def identify_taxon_hierarchy(tree, s2t):
    """
Recursively follow the taxon hierarchy of species on the tip
of the tree.

Parameters
----------
tree
    ete3 Tree instance
    with leaf names corresponding to the entries
    in species set, s2t dict and t2s dict
s2t : dict
    species indexed connections to taxa
    {key : string}
    """
    for leaf in tree.iter_leaves():
        leaf.add_feature("taxa", [])

        # each leaf should have at least one level of taxon
        # but be on the safe side and set a default value too
        leaf_regex = re.compile(leaf.name + "*")
        search_result = list(filter(leaf_regex.match, s2t.keys()))
        if search_result:
            species = list(filter(leaf_regex.match, s2t.keys()))[0]
        else:
            print("{} leaf was not found in taxaonomy file, check your input".format(leaf.name))
            exit(-1)
 
        taxon = s2t[species]

        leaf.taxa.append(taxon)

        # search for higher level taxa
        while taxon in list(s2t.keys()):
            taxon = s2t[taxon]
            leaf.taxa.append(taxon)


def check_monophyly(tree, species, logical_taxa, taxa2species, species_to_ignore):
    """
Check monophyly of the taxa described in the csv file

tree
    ete3 Tree instance
    with leaf names corresponding to the entries
    in species set, s2t dict and t2s dict
species : set
    list of species
taxa : set
    list of taxa
s2t : dict
    species indexed connections to taxa
    {key : string}
t2s : dict
    taxon indexed connections to species
    {key : [string1, string2, ..., stringN]}
    """
    for taxon in logical_taxa:
        taxon_group = []
        if [t in species for t in taxa2species[taxon]].count(False) > 0:
            # we need to collect the species recursively
            taxa_to_collect = taxa2species[taxon]

            while ([t in species for t in taxa_to_collect].count(False)) > 0:
                new_taxa_to_collect = []
                for t in taxa_to_collect:
                    if t in species:
                        # we reached species level
                        new_taxa_to_collect.append(t)
                    else:
                        # we can still expand to other taxon
                        new_taxa_to_collect.extend(
                            [s for s in taxa2species[t] if s not in species_to_ignore]
                        )
                    taxa_to_collect = new_taxa_to_collect
            taxon_group = taxa_to_collect
        else:
            # they are real species
            taxon_group = [s for s in taxa2species[taxon] if s not in species_to_ignore]

            # only consider taxa with at least two species
            if len(taxon_group) > 1:
                is_monophyletic, clade_type, violators = tree.check_monophyly(
                    values=taxon_group, target_attr="name",
                )

                for s in taxon_group:
                    leaf = tree.search_nodes(name=s)[0]
                    leaf.add_feature("clade_type", clade_type)

                if is_monophyletic:
                    print("Taxon {} is monophyletic".format(taxon))
                else:
                    # violators_in_taxon = set()
                    # print(taxa2species[taxon])
                    # for v in violators:
                    #    if v in taxa2species[taxon]:
                    #        violators_in_taxon.add(v.name)
                    print(
                        "Taxon {} is not monophyletic, but {}\n Following leaves violate it: {}".format(
                            taxon,
                            clade_type,
                            # violators_in_taxon,
                            ", ".join(sorted([v.name for v in violators])),
                        )
                    )


if __name__ == "__main__":
    main()

#https://msprime.readthedocs.io/en/stable/

import msprime

# function: msprime.simulate
# Parameters:
# sample_size (int) The number of individuals in our sample.
# Ne (float) The effective (diploid) population size
# length (float) The length of the simulated region in bases.
# recombination_rate (float) The rate of recombination per base per generation
# mutation_rate (float) The rate of mutation per base per generation
# random_seed (int) The random seed


# Simulate tree
tree_sequence = msprime.simulate(sample_size=5, Ne=1000, random_seed=10)


# Accessing the tree
tree1 = next( tree_sequence.trees() )
print('Here is the print of the tree in sparse oriented format explained by the paper:')
print(tree1)
print('')


# Drawing tree
msprime.SparseTree.draw(tree1, path = 'sampletree1.svg', width = 1024, height = 864)

# Run terminal commands directly from python in order to convert the svg format files to
# common png format (the svg format pictures provided by msprime are not very clear)
import os
os.system('convert sampletree1.svg sampletree1.png')
os.system('rm sampletree1.svg')


# Save tree in HDF or HDF5 format (zlib_compression = False or True)
msprime.TreeSequence.dump(self = tree_sequence, path = "sampletree1", zlib_compression = False)


print('Here we print the number of trees which are compacted into tree_sequence object')
print( tree_sequence.get_num_trees()   )
print('')


print('The total number of records are printed using the following command')
print( tree_sequence.get_num_records() )
print('')


print('Prints the sequence length which tree is constructed based on that. In this example = 1')
print( tree_sequence.get_sequence_length() )
print('')



print('Print the coalescence records for the tree_sequence object, including:')
for record in tree_sequence.records():
   print(record)
print('')


print('get_parent function can find the parent for a given node in the tree')
print( tree1.get_parent(1) )
print('')


print('To check what get_time function can do, lets look at the tree at sampletree1.png')
print( tree1.get_time(0) )
print( tree1.get_time(1) )
print( tree1.get_time(2) )
print( tree1.get_time(3) )
print( tree1.get_time(4) )
print( tree1.get_time(5) )
print( tree1.get_time(6) )
print( tree1.get_time(7) )
print( tree1.get_time(8) )
print('')



# We can also print the node and time in a more fashionable manner using the following code
# get_time() returns the time of the specified node in generations.
print('Node and their time: ')
u = 0
while u != msprime.NULL_NODE:
   print("node {}: time = {}".format(u, tree1.get_time(u)))
   u = tree1.get_parent(u)


print('Branch length for node 0')
print( tree1.get_branch_length(0) )
print('')


print('Get the time to the most recent common ancestor of node 0 and 1')
print( tree1.get_tmrca(0,1) )
print('')


print('Get the most recent common ancestor of node 0 and 1')
print( tree1.get_mrca(0,1) )
print('')



print('Get the time to the most recent common ancestor of node 3 and 4')
print( tree1.get_tmrca(3,4) )
print('')



# get_num_leaves() function can be used to return the number of leaves in this tree underneath
# the specified node.
print( 'Find the number of leaves underneath node 7' )
print( tree1.get_num_leaves(7) )
print('')



# get_children function returns the child of the specified node
print( 'Child of node 7 is:' )
print( tree1.get_children(7) )
print('')



# get_root function returns the root of the tree
print( 'The root of this tree is:' )
print( tree1.get_root() )
print('')



# is_internal(u)
# Returns True if the specified node is not a leaf.
# is_leaf(u)
# Returns True if the specified node is a leaf



# Recombination
# The length parameter specifies the length of the simulated sequence in bases,
# and may be a floating point number. If length is not supplied, it is assumed to be 1.
# The recombination_rate parameter specifies the rate of crossing over per base per generation,
# and is zero by default.


# Generates a tree with recombination
tree_sequence = msprime.simulate(
   sample_size=5, Ne=1000, length=1e4, recombination_rate=2e-8, random_seed=10
   )


# Since we included recombination_rate parameter, then more than one tree is generated
#  get_interval() returns a tuple (l, r) representing the left-most (inclusive) and right-most
#   (exclusive) coordinates of the genomic region covered by this tree.
print('Trees and the genomic intervals which these trees are constructed on:')
for tree in tree_sequence.trees():
   print(tree.get_interval(), str(tree) )
print('')



# Returns the number of trees inside the tree_sequence object
print('Number of trees:')
print( tree_sequence.get_num_trees() )
print('')


# Print the number of coalescent records
print('Number of records')
print( tree_sequence.get_num_records() )
print('')


# Returns the length of sequence
print('Sequence length:')
print( tree_sequence.get_sequence_length() )
print('')


# Recursively, print the coalescent records:
print('Coalescent records:')
for record in tree_sequence.records():
   print(record)
print('')


with open("records.txt", "w") as records_file:
    tree_sequence.write_records(records_file)



# Mutations
# Mutations are generated in msprime by throwing mutations down on the branches of trees at
# a particular rate. The mutations are generated under the infinite sites model, and so each
# mutation occurs at a unique (floating point) point position along the genomic 
# occupied by a tree. The mutation rate for simulations is specified using the mutation_rate
# parameter of simulate()

tree_sequence = msprime.simulate(
   sample_size=5, Ne=1000, length=1e4,
   recombination_rate=2e-8, mutation_rate=2e-8, random_seed=10)


# The get_num_mutations() function returns the number of mutations in this tree sequence.
print( "Total mutations = ", tree_sequence.get_num_mutations() )
print('')




# list function in python converts sequence types data to lists
# tree.mutations returns an object with 3 attributes
# position is the location of the mutation in genomic coordinates,
# node is the node in the tree above which the mutation occurs,
# index is the (zero-based) index of the mutation in the

print('Following is the information regarding the mutations in the tree')
for tree in tree_sequence.trees():
    print(tree.get_interval(), list(tree.mutations()))
print('')



# Write the mutations in a text format
with open("mutations.txt", "w") as mutations_file:
    tree_sequence.write_mutations(mutations_file, header = True)

print('')

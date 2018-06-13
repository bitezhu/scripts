from ete3 import Tree
import sys

t = Tree(sys.argv[1])

#print dir(t)
tnode =  t.get_children()[0]
print tnode.name
print tnode.children

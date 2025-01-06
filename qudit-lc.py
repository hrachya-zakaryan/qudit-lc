from graph import Graph

g=Graph({(1,2):2,(1,3):2,(2,3):1,(2,4):2,(3,4):1},3)
print(g.graph)
g.complement(3,1)
print(g.graph)

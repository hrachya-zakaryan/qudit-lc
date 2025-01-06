from collections import defaultdict

class Graph:
    """
    A class representing a qudit graph.
    
    Attributes:
        graph (defaultdict): The desired graph written in the form {(u,v):weight, ...}
        d (int): qudit dimension
    """

    def __init__(self, graph, d):
        """
        Initializes the class, the graph can also be passed as a dictionary and will be converted to defaultdict.
        """
        self.graph=defaultdict(int, graph) #Missing key default to weight 0
        self.d=d

    def scale(self, v, k):
        """
        Performs local scaling with scaling factor k on vertex v. weight_vu' = k*weight_vu mod d. Graph updated in place.

        Args:
            v (int): vertex number
            k (int): scaling factor, must be between 1 and d-1
        """

        if k<1 or k>self.d-1 or not isinstance(k,int):
            print(f"Please input an integer between 1 and {self.d}.")
            return
        
        vertex_passed=False #Flag for reaching v in the first position of the edge tuple

        for edge in sorted(self.graph.keys()):
            if edge[0]!=v and vertex_passed==True: #If we have passed all the edges with v in the first position there are no more edges that include v
                break
            if edge[0]==v or edge[1]==v: 
                if edge[0]==v:
                    vertex_passed=True #Turn the flag when v first seen in first position
                self.graph[edge]=self.graph[edge]*k % self.d #scale the edge connected to v
        
    def complement(self, v, k):
        """
        Performs local complementation with scaling factor k on vertex v. weight_uw' = weight_uw+k*weight_vu*weight_vw mod d. Graph updated in place.

        Args:
            v (int): vertex number
            k (int): scaling factor, must be between 1 and d-1
        """
        if k<1 or k>self.d-1 or not isinstance(k,int):
            print(f"Please input an integer between 1 and {self.d}.")
            return
        
        neighbours=self.neighbourhood(v) #List of neighbours of v

        for i in range(len(neighbours)): #Considering every neighbour pair
            for j in range(i+1, len(neighbours)):
                u=neighbours[i]
                w=neighbours[j]
                weight_uw = (self.graph[(u,w)]+k*self.graph[(min(u,v),max(u,v))]*self.graph[(min(w,v),max(w,v))]) % self.d #Calculates the new weight of the edge (u,w)
                if weight_uw==0:
                    self.graph.pop((u,w)) #Remove edge if the new weight is 0
                else:
                    self.graph[(u,w)]=weight_uw #Updates the weight, if the edge wasn't present adds the edge

    def neighbourhood(self,v):
        """
        Finds the neighbours of the vertex v.
        
        Args:
            v (int): The vertex whose neighbours are needed
        Returns:
            neighbours (list[int]): The list of neighbours of v
        """
        neighbours=[]
        for edge in sorted(self.graph.keys()):
            if edge[0]<v:
                if v in edge:
                    neighbours.append(edge[0])
                else:
                    continue
            elif edge[0]==v:
                neighbours.append(edge[1])
            else: #After all the edges where v is in the first position, no more edges include v
                break
        return neighbours